#!/usr/bin/env python3

import gzip
import json
import os
import re
import subprocess
import time

CLAIR3_CONNECT_BIN = os.environ.get("CLAIR3_CONNECT_BIN",
    "<WORKSPACE_ROOT>/target/release/clair3-connect")
SERVER = os.environ.get("CLAIR3_CONNECT_SERVER", "http://127.0.0.1:50071")
WORKDIR = os.environ.get("CLAIR3_CONNECT_AGENT_WORK", "/tmp/clair3_connect_agent")
os.makedirs(WORKDIR, exist_ok=True)

def _e(k, d=""):
    return os.environ.get(k, d)

CALLABLE_MIN_DEPTH = _e("CLAIR3_CONNECT_CALLABLE_MIN_DEPTH")
MIN_PQ = int(_e("CLAIR3_CONNECT_MIN_PQ", "20"))

_REF_BY_PLAT = {"ont": _e("REF_ONT"), "hifi": _e("REF_HIFI"), "ilmn": _e("REF_ILMN")}

SAMPLES = {}

APOE_TRUTH_JSON = _e("APOE_TRUTH_JSON",
                     os.path.join(WORKDIR, "..", "apoe", "apoe_truth.json"))

def _load_samples():
    """Load the sample registry from $CLAIR3_CONNECT_SAMPLES (JSON:
    {name: {"bam": ..., "ref": ...}}), falling back to the APOE truth file used
    by the original experiment ({name: {"bam": ..., "platform": ...}})."""
    path = _e("CLAIR3_CONNECT_SAMPLES")
    try:
        d = json.load(open(path)) if path else json.load(open(APOE_TRUTH_JSON))
    except Exception:
        return
    for s, info in d.items():
        bam = info.get("bam", "")
        if not bam:
            continue
        ref = info.get("ref") or _REF_BY_PLAT.get(info.get("platform", "ilmn"), _e("REF_ILMN"))
        SAMPLES[s] = {"bam": bam, "ref": ref}

_load_samples()

def _health(server=None):
    import sys
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "clair3_connect_pb"))
    import grpc, clair3_connect_pb2, clair3_connect_pb2_grpc
    target = (server or SERVER).replace("http://", "").replace("https://", "")
    with grpc.insecure_channel(target) as ch:
        stub = clair3_connect_pb2_grpc.Clair3ConnectInferenceStub(ch)
        return stub.Health(clair3_connect_pb2.HealthRequest(), timeout=15)

def discover_models(server: str = None) -> dict:
    """List the models the Clair3-connect server serves, with the platform and
    stage each targets as declared in the server's model registry."""
    resp = _health(server)
    models = [{"model_id": m.model_id, "platform": m.platform, "stage": m.stage,
               "description": m.description} for m in resp.models]
    source = "registry"
    if not models:
        source = "model_id_naming_convention"
        for mid in sorted(resp.registered_models.keys()):
            plat = "ont" if "ont" in mid else "hifi" if "hifi" in mid else "ilmn" if "ilmn" in mid else "?"
            stage = "pileup" if mid.startswith("pileup") else "full_alignment"
            models.append({"model_id": mid, "platform": plat, "stage": stage})
    return {"server_version": resp.server_version, "platform_source": source, "models": models}

def _fa_model_for(platform: str):
    models = discover_models()["models"]
    for m in models:
        if m["platform"] == platform and m["stage"] == "full_alignment":
            return m["model_id"], None
    avail = sorted({m["platform"] for m in models if m["stage"] == "full_alignment"})
    return None, {"error": f"no full_alignment model registered for platform '{platform}'",
                  "platforms_with_models": avail}

def _sample(sample):
    s = (sample or "").strip()
    cfg = SAMPLES.get(s)
    if not cfg:
        return s, None, {"error": f"no registered sample '{sample}'",
                         "available_samples": sorted(SAMPLES)}
    if not cfg["bam"] or not os.path.exists(cfg["bam"]):
        return s, None, {"error": f"BAM not available for sample '{sample}'"}
    return s, cfg, None

def bam_header(sample: str) -> dict:
    """Return a registered sample's BAM header records — the @HD/@RG/@PG/@CO lines
    `samtools view -H` shows (the long @SQ contig list is summarized to a count),
    so the agent can DETERMINE the sequencing platform ITSELF from the evidence
    (e.g. @RG `PL:`, the aligner named in @PG, read-group/flowcell hints) and then
    pick the matching model. The platform is deliberately NOT labeled here. The
    sample is addressed by name; the BAM path stays server-side. `warnings`
    states header facts that conflict with the registration (an @RG SM tag other
    than the registered name) or that leave the platform without evidence."""
    s, cfg, err = _sample(sample)
    if err:
        return err
    samtools = os.environ.get("SAMTOOLS",
        "<CONDA>/envs/clair3-torch/bin/samtools")
    try:
        p = subprocess.run([samtools, "view", "--no-PG", "-H", cfg["bam"]],
                           capture_output=True, text=True, timeout=120)
    except Exception as e:
        return {"error": f"samtools failed: {e}"}
    if p.returncode != 0:
        return {"error": "samtools view -H failed", "stderr_tail": p.stderr[-300:]}
    lines = p.stdout.splitlines()
    keep = [ln[:220] for ln in lines if ln.startswith(("@HD", "@RG", "@PG", "@CO"))][:60]
    n_sq = sum(1 for ln in lines if ln.startswith("@SQ"))
    rg = [dict(f.split(":", 1) for f in ln.split("\t")[1:] if ":" in f)
          for ln in lines if ln.startswith("@RG")]
    rg_samples = sorted({r["SM"] for r in rg if r.get("SM")})
    rg_platforms = sorted({r["PL"] for r in rg if r.get("PL")})
    warnings = []
    if rg_samples and s not in rg_samples:
        warnings.append(f"@RG SM tag(s) {rg_samples} differ from the registered sample name '{s}'")
    if not rg_platforms and not any(ln.startswith("@PG") for ln in lines):
        warnings.append("header has no @RG PL tag and no @PG records: it carries no evidence "
                        "of the sequencing platform")
    res = {"sample": s, "header": keep, "n_reference_contigs": n_sq,
           "read_group_samples": rg_samples, "read_group_platforms": rg_platforms,
           "note": "infer the sequencing platform from these records, then select "
                   "the matching model (see discover_models)"}
    if warnings:
        res["warnings"] = warnings
    return res

_REGION_RE = re.compile(r"^[^:\s]+:\d+-\d+$")

def call_variants(sample: str, platform: str, region: str, threads: int = 8) -> dict:
    """Call germline variants for a registered sample over a region with the
    full-alignment model registered for `platform`. Returns the output VCF path,
    a summary of PASS variants, and the region's callability (share of bases
    with enough well-mapped coverage for Clair3 to assess them)."""
    s, cfg, err = _sample(sample)
    if err:
        return err
    platform = (platform or "").strip().lower()
    region = (region or "").strip().replace(",", "")
    if not _REGION_RE.match(region):
        return {"error": f"region '{region}' is not chr:start-end"}
    model, err = _fa_model_for(platform)
    if err:
        return err
    out = os.path.join(WORKDIR, f"{s}_{platform}_{region.replace(':', '_')}")
    os.makedirs(out, exist_ok=True)
    cmd = ["env", "-u", "LD_LIBRARY_PATH", CLAIR3_CONNECT_BIN,
           "--bam", cfg["bam"], "--ref", cfg["ref"], "--region", region,
           "--fa-model", model, "--server", SERVER, "--threads", str(threads),
           "--platform", platform, "--sample-name", s, "--output", out, "--non-interactive",
           "--callability"]
    if CALLABLE_MIN_DEPTH:
        cmd += ["--callable-min-depth", CALLABLE_MIN_DEPTH]
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
    res = {"vcf_path": vcf, "region": region, "sample": s, "platform": platform,
           "model_id": model, "elapsed_s": round(dt, 2), **summ}
    try:
        cj = json.load(open(os.path.join(out, "callability.json")))
        res["callability"] = {"callable_fraction": cj["callable_fraction"],
                              "median_depth": cj["median_depth"],
                              "bp_by_state": cj["bp_by_state"],
                              "min_depth": cj["params"]["min_depth"]}
        res["callable_bed"] = os.path.join(out, "callable.bed")
    except Exception:
        pass
    return res

def _open(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)

def _summarize(vcf_path: str, region: str = None) -> dict:
    """Count SNV/indel and het/hom among PASS non-reference records."""
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
            gt = dict(zip(c[8].split(":"), c[9].split(":"))).get("GT", "./.").replace("|", "/")
            alleles = gt.split("/")
            if c[6] not in ("PASS", ".") or not any(a not in ("0", ".") for a in alleles):
                continue
            rc += 1
            is_snv = len(ref) == 1 and all(len(a) == 1 for a in alt.split(","))
            rs += 1 if is_snv else 0
            ri += 0 if is_snv else 1
            if len(alleles) == 2 and alleles[0] == alleles[1]:
                hom += 1
            else:
                het += 1
    return {"total": rc, "snv": rs, "indel": ri, "het": het, "hom": hom,
            "het_fraction": round(het / rc, 4) if rc else None}

def _run_json(cmd, timeout=600):
    p = subprocess.run(["env", "-u", "LD_LIBRARY_PATH"] + cmd,
                       capture_output=True, text=True, timeout=timeout)
    if p.returncode != 0:
        return None, {"error": f"{cmd[1]} failed", "stderr_tail": p.stderr[-600:]}
    try:
        return json.loads(p.stdout), None
    except json.JSONDecodeError:
        return None, {"error": f"{cmd[1]} returned unreadable output"}

def _vcf_records(vcf_path):
    recs = []
    with _open(vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 10:
                continue
            gt = dict(zip(c[8].split(":"), c[9].split(":"))).get("GT", "./.").replace("|", "/")
            recs.append({"chrom": c[0], "pos": int(c[1]), "ref": c[3], "alt": c[4],
                         "filter": c[6], "gt": gt, "qual": c[5]})
    return recs

def _locus_state(bed_path, chrom, pos):
    with open(bed_path) as f:
        for line in f:
            c = line.split("\t")
            if c[0] == chrom and int(c[1]) < pos <= int(c[2]):
                return c[3].strip()
    return None

def genotypes_at(sample: str, platform: str, region: str, positions: list,
                 include_evidence: bool = False) -> dict:
    """Call variants for (sample, platform) over `region` with Clair3, then
    report each requested position as one of three states:
      variant_called   - Clair3 emitted a PASS variant record (its GT is given);
      callable_hom_ref - Clair3 emitted no variant and the locus was CALLABLE
                         (enough well-mapped reads for Clair3 to assess it);
      no_call          - the genotype is unknown; `reasons` says why
                         (region_not_processed, no_coverage, low_coverage,
                         poor_mapping_quality, excessive_coverage,
                         low_quality_call).
    Genotypes come from Clair3 only. With include_evidence=true each position
    also carries descriptive read evidence (depth, MAPQ, base counts).
    `positions` are 'chrom:pos' strings (1-based)."""
    if not isinstance(positions, list) or not positions:
        return {"error": "positions must be a non-empty list of 'chrom:pos'"}
    r = call_variants(sample, platform, region)
    if "error" in r:
        return r
    recs = _vcf_records(r["vcf_path"])
    rchr, se = r["region"].split(":")
    rs, re_ = (int(x) for x in se.split("-"))
    out = []
    for p in positions:
        try:
            chrom, pos = p.strip().split(":")[:2]
            pos = int(pos.replace(",", ""))
        except (ValueError, AttributeError):
            out.append({"query": p, "error": "expected 'chrom:pos'"})
            continue
        item = {"position": f"{chrom}:{pos}"}
        exact = next((x for x in recs if x["chrom"] == chrom and x["pos"] == pos), None)
        spanning = next((x for x in recs if x["chrom"] == chrom and x["pos"] < pos < x["pos"] + len(x["ref"])
                         and x["filter"] in ("PASS", ".") and x["gt"] not in ("0/0", "./.")), None)
        non_ref = lambda x: any(a not in ("0", ".") for a in x["gt"].split("/"))
        locus = _locus_state(r["callable_bed"], chrom, pos) if r.get("callable_bed") else None
        if chrom != rchr or not rs <= pos <= re_:
            item.update(state="no_call", genotype=None, reasons=["region_not_processed"])
        elif exact and non_ref(exact) and exact["filter"] in ("PASS", "."):
            item.update(state="variant_called", genotype=exact["gt"], ref=exact["ref"], alt=exact["alt"],
                        reasons=[])
        elif spanning:
            item.update(state="variant_called", genotype=spanning["gt"], ref=spanning["ref"],
                        alt=spanning["alt"], reasons=[f"inside_called_deletion_at_{spanning['pos']}"])
        elif exact and non_ref(exact):
            item.update(state="no_call", genotype=None, reasons=["low_quality_call"],
                        clair3_record={k: exact[k] for k in ("ref", "alt", "filter", "gt", "qual")})
        elif locus == "CALLABLE":
            item.update(state="callable_hom_ref", genotype="0/0", reasons=[])
        else:
            item.update(state="no_call", genotype=None,
                        reasons=[(locus or "region_not_processed").lower()])
        if locus:
            item["locus_state"] = locus
        out.append(item)
    res = {"sample": r["sample"], "platform": r["platform"], "model_id": r["model_id"],
           "region": r["region"], "vcf_path": r["vcf_path"],
           "callability": r.get("callability"), "genotypes": out,
           "note": "genotypes are Clair3's; no_call means unknown, not homozygous reference; "
                   "genotypes are unphased"}
    if include_evidence:
        cfg = SAMPLES[r["sample"]]
        ok = [g["position"] for g in out if "position" in g]
        d, err = _run_json([CLAIR3_CONNECT_BIN, "sites", "--bam", cfg["bam"], "--ref", cfg["ref"],
                            "--sites", ",".join(ok)])
        if err:
            res["evidence_error"] = err
        else:
            ev = {f"{x['contig']}:{x['pos']}": x["evidence"] for x in d["sites"]}
            for g in out:
                if g.get("position") in ev:
                    g["evidence"] = ev[g["position"]]
    return res

def phase_at(sample: str, positions: list) -> dict:
    """Read-backed phasing of heterozygous SNVs. Counts the reads/fragments that
    span each pair of positions and reports whether the ALT alleles lie on the
    same haplotype (cis) or on different haplotypes (trans), with a phred-scaled
    phase quality (PQ). status 'phased' means every position is in one block of
    confident relations; otherwise the reads do not determine the phase.
    `positions` are 'chrom:pos:REF:ALT' (ALT is taken from the reads if only
    'chrom:pos' is given)."""
    if not isinstance(positions, list) or len(positions) < 2:
        return {"error": "positions must list at least two 'chrom:pos:REF:ALT' sites"}
    s, cfg, err = _sample(sample)
    if err:
        return err
    d, err = _run_json([CLAIR3_CONNECT_BIN, "phase", "--bam", cfg["bam"], "--ref", cfg["ref"],
                        "--sites", ",".join(p.strip() for p in positions),
                        "--min-pq", str(MIN_PQ)])
    if err:
        return err
    d["sample"] = s
    return d

PHARMVAR_FILES = {}
_REFSEQ_CHROM = {"NC_000010.11": "chr10", "NC_000022.11": "chr22", "NC_000001.11": "chr1",
                 "NC_000007.14": "chr7", "NC_000006.12": "chr6", "NC_000019.10": "chr19"}
STRUCTURAL_CAVEATS = {
    "CYP2D6": "CYP2D6 star alleles also include a whole-gene deletion (*5), gene duplications/multiplications "
              "(xN) and CYP2D6/CYP2D7 hybrid genes (e.g. *13, *36, *68). These are copy-number/structural "
              "changes that are not in the small-variant definitions below and cannot be assessed from small-"
              "variant calls; a diplotype from small variants alone is not reliable for CYP2D6 without a "
              "copy-number/structural analysis (for this reason PharmCAT, for example, does not call CYP2D6 "
              "from a VCF).",
    "CYP2A6": "CYP2A6 includes deletion, duplication and CYP2A6/CYP2A7 hybrid alleles not assessable from "
              "small-variant calls.",
}

def _pharmvar_file(gene):
    if gene in PHARMVAR_FILES:
        return PHARMVAR_FILES[gene]
    import glob
    root = _e("CLAIR3_CONNECT_PHARMVAR_DIR")
    hits = sorted(glob.glob(os.path.join(root, f"{gene}-*", "GRCh38", "*.haplotypes.tsv"))) if root else []
    return hits[-1] if hits else None

def allele_definitions(gene: str, alleles: list = None, max_alleles: int = 60) -> dict:
    """PharmVar star-allele definitions for a pharmacogene, as structured data:
    the core allele definitions (sub-alleles collapsed onto their core allele),
    which allele the GRCh38 reference sequence itself carries, every defining
    position and one region covering them (for genotypes_at), and gene-level
    caveats (alleles a small-variant caller cannot assess). Definitions only:
    genotypes come from genotypes_at, and matching them to alleles is up to you."""
    gene = (gene or "").strip().upper()
    path = _pharmvar_file(gene)
    if not path or not os.path.exists(path):
        return {"error": f"no PharmVar definitions available for '{gene}'"}
    version, rows = None, []
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                version = line.strip().lstrip("#")
                continue
            c = line.rstrip("\n").split("\t")
            if c[0] == "Haplotype Name" or len(c) < 9:
                continue
            rows.append(c)
    core_rows = {}
    subs = {}
    reference_allele = None
    for c in rows:
        name = c[0]
        core = name.split(".")[0]
        star = "*" + core.split("*", 1)[1]
        subs.setdefault(star, set())
        if name != core:
            subs[star].add(name)
            continue
        if c[3] == "REFERENCE":
            reference_allele = star
            core_rows.setdefault(star, [])
            continue
        chrom = _REFSEQ_CHROM.get(c[3], c[3])
        start, ref, alt, vtype = int(c[4]), c[6] or "-", c[7] or "-", c[8] or "?"
        core_rows.setdefault(star, []).append(
            f"{c[2] or '-'} {chrom}:{start} {ref}>{alt}" + ("" if vtype == "substitution" else f" ({vtype})"))
    order = sorted(core_rows, key=lambda a: (int("".join(ch for ch in a[1:] if ch.isdigit()) or 0), a))
    want = [("*" + a.split("*")[-1]).upper() for a in alleles] if alleles else None
    shown = [a for a in order if (want is None or a.upper() in want)]
    positions = sorted({v.split()[1] for a in core_rows for v in core_rows[a]},
                       key=lambda p: (p.split(":")[0], int(p.split(":")[1])))
    region = None
    if positions:
        ch = positions[0].split(":")[0]
        ps = [int(p.split(":")[1]) for p in positions if p.startswith(ch + ":")]
        region = f"{ch}:{max(1, min(ps) - 100)}-{max(ps) + 100}"
    res = {"gene": gene, "source": f"PharmVar {version or ''} GRCh38 haplotype definitions".strip(),
           "caveats": [STRUCTURAL_CAVEATS[gene]] if gene in STRUCTURAL_CAVEATS else [],
           "reference_allele": (f"the GRCh38 reference sequence carries {gene}{reference_allele}: a sample with no "
                                f"variant at any defining position is {reference_allele}/{reference_allele}, and "
                                f"alleles are defined relative to it") if reference_allele else None,
           "n_core_alleles": len(order),
           "n_suballeles_collapsed": sum(len(v) for v in subs.values()),
           "region_covering_positions": region,
           "defining_positions": positions if len(positions) <= 120 else f"{len(positions)} positions",
           "variant_format": "rsID chrom:pos REF>ALT (PharmVar coordinates; '-' = deleted/inserted bases, "
                             "an indel appears in a VCF at pos-1 with an anchor base)",
           "note": "definitions only; get genotypes with genotypes_at (Clair3) and phase with phase_at when "
                   "an allele needs variants on the same haplotype",
           "core_alleles": {}}
    for a in shown[:max_alleles]:
        res["core_alleles"][a] = core_rows[a]
        if len(json.dumps(res)) > 5000:
            del res["core_alleles"][a]
            break
    left = len(shown) - len(res["core_alleles"])
    if left > 0:
        res["truncated"] = f"{left} more core alleles not shown; pass `alleles` to request specific ones"
    return res

DISPATCH = {
    "discover_models": discover_models,
    "bam_header": bam_header,
    "call_variants": call_variants,
    "genotypes_at": genotypes_at,
    "phase_at": phase_at,
}

TOOLS = [
    {"type": "function", "function": {
        "name": "discover_models",
        "description": "List variant-calling models the Clair3-connect server serves, with the platform (ont/hifi/ilmn) and stage each targets as declared in the server's model registry.",
        "parameters": {"type": "object", "properties": {}, "required": []}}},
    {"type": "function", "function": {
        "name": "bam_header",
        "description": "Return a registered sample's BAM header records (@HD/@RG/@PG/@CO, as `samtools view -H` shows; @SQ contigs summarized to a count). Use it to determine the sample's sequencing platform yourself (from @RG PL, the aligner in @PG, flowcell/read-group hints) and then pick the matching model. Platform is not pre-labeled. Also lists the @RG sample (SM) and platform (PL) tags; `warnings` flags an SM tag that differs from the registered sample name, or a header with no platform evidence.",
        "parameters": {"type": "object", "properties": {
            "sample": {"type": "string", "description": "registered sample name, e.g. HG00136"}},
            "required": ["sample"]}}},
    {"type": "function", "function": {
        "name": "call_variants",
        "description": "Call germline SNVs/indels for a registered sample over a genomic region using the Clair3-connect full-alignment model for the given platform. Returns the VCF path, counts of PASS variants, and the region's callability (share of bases with enough well-mapped coverage for Clair3 to assess them, median depth).",
        "parameters": {"type": "object", "properties": {
            "sample": {"type": "string", "description": "e.g. HG00136"},
            "platform": {"type": "string", "enum": ["ont", "hifi", "ilmn"]},
            "region": {"type": "string", "description": "chr:start-end, e.g. chr19:44908000-44909000"},
        }, "required": ["sample", "platform", "region"]}}},
    {"type": "function", "function": {
        "name": "genotypes_at",
        "description": "Call variants for (sample, platform) over a region with Clair3, then report each requested 1-based position as one of three states: variant_called (Clair3 PASS variant, with its GT), callable_hom_ref (no Clair3 variant at a CALLABLE locus, i.e. enough well-mapped reads for Clair3 to assess it), or no_call (genotype unknown; reasons such as region_not_processed, low_coverage, no_coverage, poor_mapping_quality, excessive_coverage, low_quality_call). A no_call is NOT a reference genotype. Genotypes are unphased. Set include_evidence=true to also get descriptive read evidence per position (depth, MAPQ, base counts).",
        "parameters": {"type": "object", "properties": {
            "sample": {"type": "string"}, "platform": {"type": "string", "enum": ["ont", "hifi", "ilmn"]},
            "region": {"type": "string", "description": "chr:start-end covering the positions"},
            "positions": {"type": "array", "items": {"type": "string"},
                          "description": "list of 'chrom:pos', e.g. ['chr19:44908684','chr19:44908822']"},
            "include_evidence": {"type": "boolean", "description": "also return read evidence (depth, MAPQ, base counts) per position; default false"},
        }, "required": ["sample", "platform", "region", "positions"]}}},
    {"type": "function", "function": {
        "name": "phase_at",
        "description": "Read-backed phasing of heterozygous SNVs from reads/fragments spanning two or more positions. Reports per-pair cis/trans counts with a phred-scaled phase quality (PQ) and the resulting haplotypes; status 'phased' means all positions are in one confident block, otherwise the phase is not determined by the reads.",
        "parameters": {"type": "object", "properties": {
            "sample": {"type": "string"},
            "positions": {"type": "array", "items": {"type": "string"},
                          "description": "list of 'chrom:pos:REF:ALT' SNVs, e.g. ['chr19:44908684:T:C','chr19:44908822:C:T']"},
        }, "required": ["sample", "positions"]}}},
]

TOOLS_CLINICAL = TOOLS

DISPATCH_PGX = {"allele_definitions": allele_definitions}
PGX_TOOLS = [
    {"type": "function", "function": {
        "name": "allele_definitions",
        "description": "PharmVar star-allele definitions for a pharmacogene as structured data: core allele definitions (sub-alleles collapsed), which allele the GRCh38 reference itself carries, all defining positions and a region covering them (to pass to genotypes_at), and gene-level caveats such as alleles a small-variant caller cannot assess. Definitions only: genotypes come from genotypes_at; matching genotypes to alleles is up to you.",
        "parameters": {"type": "object", "properties": {
            "gene": {"type": "string", "description": "e.g. CYP2C19"},
            "alleles": {"type": "array", "items": {"type": "string"},
                        "description": "optional: only these core alleles, e.g. ['*2','*17']"},
        }, "required": ["gene"]}}},
]

if __name__ == "__main__":
    import sys
    print(json.dumps([t["function"]["name"] for t in TOOLS]))
    if len(sys.argv) > 1 and sys.argv[1] == "discover":
        print(json.dumps(discover_models(), indent=2))
