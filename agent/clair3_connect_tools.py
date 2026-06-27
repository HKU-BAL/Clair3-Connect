#!/usr/bin/env python3

import gzip
import json
import os
import subprocess
import time

CLAIR3_CONNECT_BIN = os.environ.get("CLAIR3_CONNECT_BIN",
    "<WORKSPACE_ROOT>/target/release/clair3-connect")
SERVER = os.environ.get("CLAIR3_CONNECT_SERVER", "http://127.0.0.1:50071")
WORKDIR = os.environ.get("CLAIR3_CONNECT_AGENT_WORK", "/tmp/clair3_connect_agent")
os.makedirs(WORKDIR, exist_ok=True)

def _e(k, d=""):
    return os.environ.get(k, d)

_MODELS_BY_PLAT = {
    "ont":  ("pileup_ont_hac_v500", "full_alignment_ont_hac_v500"),
    "hifi": ("pileup_hifi_revio",   "full_alignment_hifi_revio"),
    "ilmn": ("pileup_ilmn",         "full_alignment_ilmn"),
}
_REF_BY_PLAT = {"ont": _e("REF_ONT"), "hifi": _e("REF_HIFI"), "ilmn": _e("REF_ILMN")}

SAMPLES = {}

APOE_TRUTH_JSON = _e("APOE_TRUTH_JSON",
                     os.path.join(WORKDIR, "..", "apoe", "apoe_truth.json"))

def _load_apoe_samples():
    try:
        d = json.load(open(APOE_TRUTH_JSON))
    except Exception:
        return
    for s, info in d.items():
        plat = info.get("platform", "ilmn")
        bam = info.get("bam", "")
        if not bam:
            continue
        pm, fm = _MODELS_BY_PLAT.get(plat, _MODELS_BY_PLAT["ilmn"])
        SAMPLES[(s, plat)] = {"bam": bam, "ref": _REF_BY_PLAT.get(plat, _e("REF_ILMN")),
                              "pileup_model": pm, "fa_model": fm}

_load_apoe_samples()

def discover_models(server: str = None) -> dict:
    """List the model_ids the Clair3-connect server currently serves (Health RPC),
    with the platform each targets. No file paths involved."""
    server = server or SERVER
    import sys
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "clair3_connect_pb"))
    import grpc, clair3_connect_pb2, clair3_connect_pb2_grpc
    target = server.replace("http://", "").replace("https://", "")
    with grpc.insecure_channel(target) as ch:
        stub = clair3_connect_pb2_grpc.Clair3ConnectInferenceStub(ch)
        resp = stub.Health(clair3_connect_pb2.HealthRequest(), timeout=15)
    models = []
    for mid in sorted(resp.registered_models.keys()):
        plat = "ont" if "ont" in mid else "hifi" if "hifi" in mid else "ilmn" if "ilmn" in mid else "?"
        stage = "pileup" if mid.startswith("pileup") else "full_alignment"
        models.append({"model_id": mid, "platform": plat, "stage": stage})
    return {"server_version": resp.server_version, "models": models}

def bam_header(sample: str) -> dict:
    """Return a registered sample's BAM header records — the @HD/@RG/@PG/@CO lines
    `samtools view -H` shows (the long @SQ contig list is summarized to a count),
    so the agent can DETERMINE the sequencing platform ITSELF from the evidence
    (e.g. @RG `PL:`, the aligner named in @PG, read-group/flowcell hints) and then
    pick the matching model. The platform is deliberately NOT labeled here. The
    sample is addressed by name; the BAM path stays server-side."""
    s = sample.strip()
    hits = [cfg for (samp, _plat), cfg in SAMPLES.items() if samp == s]
    if not hits:
        return {"error": f"no registered sample '{sample}'",
                "available_samples": sorted({samp for (samp, _p) in SAMPLES})}
    bam = hits[0].get("bam", "")
    if not bam or not os.path.exists(bam):
        return {"error": f"BAM not available for sample '{sample}'"}
    samtools = os.environ.get("SAMTOOLS",
        "<CONDA>/envs/clair3-torch/bin/samtools")
    try:
        p = subprocess.run([samtools, "view", "-H", bam],
                           capture_output=True, text=True, timeout=120)
    except Exception as e:
        return {"error": f"samtools failed: {e}"}
    if p.returncode != 0:
        return {"error": "samtools view -H failed", "stderr_tail": p.stderr[-300:]}
    lines = p.stdout.splitlines()
    keep = [ln[:220] for ln in lines if ln.startswith(("@HD", "@RG", "@PG", "@CO"))][:60]
    n_sq = sum(1 for ln in lines if ln.startswith("@SQ"))
    return {"sample": s, "header": keep, "n_reference_contigs": n_sq,
            "note": "infer the sequencing platform from these records, then select "
                    "the matching model (see discover_models)"}

def call_variants(sample: str, platform: str, region: str, threads: int = 8) -> dict:
    """Call germline variants for a registered (sample, platform) over a region
    via the Clair3-connect client/server. Returns the output VCF path + a summary."""
    key = (sample.strip(), platform.strip().lower())
    if key not in SAMPLES:
        return {"error": f"no dataset for {key}", "available": [list(k) for k in SAMPLES]}
    cfg = SAMPLES[key]
    if not cfg["bam"] or not os.path.exists(cfg["bam"]):
        return {"error": f"BAM not available for {key}"}
    out = os.path.join(WORKDIR, f"{sample}_{platform}_{region.replace(':','_')}")
    os.makedirs(out, exist_ok=True)
    cmd = ["env", "-u", "LD_LIBRARY_PATH", CLAIR3_CONNECT_BIN,
           "--bam", cfg["bam"], "--ref", cfg["ref"], "--region", region,
           "--pileup-model", cfg["pileup_model"], "--fa-model", cfg["fa_model"],
           "--server", SERVER, "--threads", str(threads), "--platform", platform,
           "--output", out, "--non-interactive"]
    t0 = time.time()
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
    dt = time.time() - t0
    vcf = os.path.join(out, "merge_output.vcf")
    if not os.path.exists(vcf):
        vcf += ".gz"
    if p.returncode != 0 or not os.path.exists(vcf):
        return {"error": "variant calling failed", "returncode": p.returncode,
                "stderr_tail": p.stderr[-800:]}
    summ = _summarize(vcf)
    return {"vcf_path": vcf, "region": region, "sample": sample, "platform": platform,
            "elapsed_s": round(dt, 2), **summ}

def _open(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)

def _summarize(vcf_path: str, region: str = None) -> dict:
    """Count SNV/indel and het/hom in a VCF (optionally a sub-region)."""
    rc = rs = ri = het = hom = 0
    rchr = rstart = rend = None
    if region:
        rchr, se = region.split(":"); rstart, rend = (int(x) for x in se.split("-"))
    with _open(vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 10:
                continue
            chrom, pos, ref, alt = c[0], int(c[1]), c[3], c[4]
            if region and (chrom != rchr or pos < rstart or pos > rend):
                continue
            rc += 1
            is_snv = len(ref) == 1 and all(len(a) == 1 for a in alt.split(","))
            rs += 1 if is_snv else 0
            ri += 0 if is_snv else 1
            gt = dict(zip(c[8].split(":"), c[9].split(":"))).get("GT", "./.").replace("|", "/")
            alleles = gt.split("/")
            if len(alleles) == 2 and alleles[0] != "." and alleles[1] != ".":
                if alleles[0] == alleles[1]:
                    hom += 1
                else:
                    het += 1
    return {"total": rc, "snv": rs, "indel": ri, "het": het, "hom": hom,
            "het_fraction": round(het / rc, 4) if rc else None}

def genotypes_at(sample: str, platform: str, region: str, positions: list) -> dict:
    """Call variants for (sample, platform) over `region`, then report the
    genotype (GT) at each requested position. A position with NO variant in the
    VCF is homozygous-reference (GT 0/0) — reported as such, not omitted, so the
    caller can reason about reference alleles. `positions` is a list of
    'chrom:pos' strings (1-based). Use to read out the two APOE-defining SNPs and
    decide a diplotype yourself."""
    r = call_variants(sample, platform, region)
    if "error" in r:
        return r
    calls = {}
    with _open(r["vcf_path"]) as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 10:
                continue
            gt = dict(zip(c[8].split(":"), c[9].split(":"))).get("GT", "./.").replace("|", "/")
            calls[(c[0], int(c[1]))] = {"ref": c[3], "alt": c[4], "gt": gt}
    out = []
    for p in positions:
        try:
            chrom, pos = p.split(":")
            pos = int(pos)
        except (ValueError, AttributeError):
            out.append({"query": p, "error": "expected 'chrom:pos'"})
            continue
        hit = calls.get((chrom, pos))
        if hit:
            out.append({"position": p, "ref": hit["ref"], "alt": hit["alt"],
                        "gt": hit["gt"], "called": True})
        else:
            out.append({"position": p, "alt": ".", "gt": "0/0", "called": False,
                        "note": "no variant called at this position — homozygous reference"})
    return {"sample": sample, "platform": platform, "region": region,
            "vcf_path": r["vcf_path"], "genotypes": out}

DISPATCH = {
    "discover_models": discover_models,
    "bam_header": bam_header,
    "call_variants": call_variants,
    "genotypes_at": genotypes_at,
}

TOOLS = [
    {"type": "function", "function": {
        "name": "discover_models",
        "description": "List variant-calling models the Clair3-connect server serves, and the platform (ont/hifi/ilmn) each targets.",
        "parameters": {"type": "object", "properties": {}, "required": []}}},
    {"type": "function", "function": {
        "name": "bam_header",
        "description": "Return a registered sample's BAM header records (@HD/@RG/@PG/@CO, as `samtools view -H` shows; @SQ contigs summarized to a count). Use it to determine the sample's sequencing platform yourself (from @RG PL, the aligner in @PG, flowcell/read-group hints) and then pick the matching model. Platform is not pre-labeled.",
        "parameters": {"type": "object", "properties": {
            "sample": {"type": "string", "description": "registered sample name, e.g. HG00136"}},
            "required": ["sample"]}}},
    {"type": "function", "function": {
        "name": "call_variants",
        "description": "Call germline SNVs/indels for a registered (sample, platform) over a genomic region using Clair3-connect. Returns the VCF path and counts.",
        "parameters": {"type": "object", "properties": {
            "sample": {"type": "string", "description": "e.g. HG00136"},
            "platform": {"type": "string", "enum": ["ont", "hifi", "ilmn"]},
            "region": {"type": "string", "description": "chr:start-end, e.g. chr19:44908000-44909000"},
        }, "required": ["sample", "platform", "region"]}}},
    {"type": "function", "function": {
        "name": "genotypes_at",
        "description": "Call variants for (sample, platform) over a region, then read out the genotype (GT) at specific 1-based positions. Positions absent from the VCF are reported as homozygous-reference (0/0). Use to inspect the APOE SNP sites and decide a diplotype.",
        "parameters": {"type": "object", "properties": {
            "sample": {"type": "string"}, "platform": {"type": "string", "enum": ["ont", "hifi", "ilmn"]},
            "region": {"type": "string", "description": "chr:start-end covering the positions"},
            "positions": {"type": "array", "items": {"type": "string"},
                          "description": "list of 'chrom:pos', e.g. ['chr19:44908684','chr19:44908822']"},
        }, "required": ["sample", "platform", "region", "positions"]}}},
]

TOOLS_CLINICAL = TOOLS

if __name__ == "__main__":
    import sys
    print(json.dumps([t["function"]["name"] for t in TOOLS]))
    if len(sys.argv) > 1 and sys.argv[1] == "discover":
        print(json.dumps(discover_models(), indent=2))
