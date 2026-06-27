#!/usr/bin/env python
"""convert_and_register.py — one command: Clair3 model -> server-ready model.

Takes a Clair3 PyTorch model directory (a dir containing `pileup.pt` and/or
`full_alignment.pt`), exports each stage to a frozen TorchScript `.ts`, and
upserts a matching entry into a server model-registry YAML — so the model is
converted, stored, and immediately servable in a single step.

Prerequisites (the conversion itself needs PyTorch + the Clair3 model code):
  * Run inside a PyTorch env (e.g. the `clair3-torch` conda env).
  * A Clair3 source clone, located via --clair3 or $CLAIR3_HOME.

Example:
  python tools/convert_and_register.py \\
      --model-dir /data/clair3_models/r1041_e82_hac_v500 \\
      --platform ont --name ont_hac_v500 --accuracy hac \\
      --registry models/registry.yaml --models-dir models \\
      --clair3 /path/to/Clair3

Produces (for --name ont_hac_v500, --stage both):
  models/pileup_ont_hac_v500.ts
  models/full_alignment_ont_hac_v500.ts
and adds/updates `pileup_ont_hac_v500` + `full_alignment_ont_hac_v500` in the
registry YAML (idempotent on model_id). Manage entries later with the server's
interactive menu: `clair3-connect-server --interactive`.

If you ALREADY have a `.ts` file, you don't need this — just register it with the
server's interactive model menu (it auto-computes the sha256).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

try:
    import yaml
except ImportError:
    sys.exit("[FATAL] PyYAML is required: `pip install pyyaml` (or use the clair3-torch env).")

HERE = Path(__file__).resolve().parent
EXPORTER = HERE / "export_torchscript" / "export_clair3_torchscript.py"
STAGES = ("pileup", "full_alignment")

def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

def build_entry(meta: dict, model_id: str, platform: str, stage: str,
                accuracy: str | None, ts_path_rel: str, ts_sha: str) -> dict:
    """Assemble a registry ModelEntry from the exporter's meta.json sidecar."""
    cfg = meta.get("config", {})
    inp = meta["input"]
    out = meta["output"]
    entry: dict = {
        "model_id": model_id,
        "platform": platform,
        "stage": stage,
    }
    if accuracy:
        entry["accuracy"] = accuracy
    entry["torchscript"] = {"path": ts_path_rel, "sha256": ts_sha}
    entry["input"] = {"shape": list(inp["shape"]), "dtype": inp["dtype"]}
    if inp.get("layout"):
        entry["input"]["layout"] = inp["layout"]
    entry["output"] = {"shape": list(out["shape"]), "dtype": out["dtype"]}
    flags = {k: cfg[k] for k in ("add_indel_length", "predict", "enable_dwell_time") if k in cfg}
    if flags:
        entry["flags"] = flags
    return entry

def upsert(registry: dict, entry: dict) -> str:
    """Add or replace an entry by model_id. Returns 'added' or 'updated'."""
    models = registry.setdefault("models", [])
    for i, m in enumerate(models):
        if m.get("model_id") == entry["model_id"]:
            models[i] = entry
            return "updated"
    models.append(entry)
    return "added"

def main() -> int:
    ap = argparse.ArgumentParser(description="Convert a Clair3 model to TorchScript and register it.")
    ap.add_argument("--model-dir", required=True,
                    help="Clair3 PyTorch model dir (contains pileup.pt and/or full_alignment.pt).")
    ap.add_argument("--platform", required=True, choices=["ont", "hifi", "ilmn"])
    ap.add_argument("--name", required=True,
                    help="Suffix for the model_id, e.g. ont_hac_v500 / hifi_revio / ilmn "
                         "(produces pileup_<name> + full_alignment_<name>).")
    ap.add_argument("--registry", required=True,
                    help="Registry YAML to upsert into (created if missing).")
    ap.add_argument("--models-dir", default=None,
                    help="Where the .ts files land (default: the registry's directory).")
    ap.add_argument("--clair3", default=None,
                    help="Path to the Clair3 source clone (else uses $CLAIR3_HOME).")
    ap.add_argument("--accuracy", default=None, help="Optional accuracy tag, e.g. hac.")
    ap.add_argument("--stage", choices=["pileup", "full_alignment", "both"], default="both")
    ap.add_argument("--trace-batch", type=int, default=4)
    args = ap.parse_args()

    if not EXPORTER.is_file():
        sys.exit(f"[FATAL] exporter not found at {EXPORTER}")
    if args.clair3:
        os.environ["CLAIR3_HOME"] = str(Path(args.clair3).resolve())

    registry_path = Path(args.registry).resolve()
    registry_dir = registry_path.parent
    models_dir = Path(args.models_dir).resolve() if args.models_dir else registry_dir
    models_dir.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="cc_convert_") as tmp:
        cmd = [sys.executable, str(EXPORTER),
               "--model-dir", str(Path(args.model_dir).resolve()),
               "--platform", args.platform,
               "--output-dir", tmp,
               "--trace-batch", str(args.trace_batch)]
        if args.stage == "pileup":
            cmd += ["--skip", "full_alignment"]
        elif args.stage == "full_alignment":
            cmd += ["--skip", "pileup"]
        print(f"[convert] running exporter: {' '.join(cmd)}")
        rc = subprocess.call(cmd)
        if rc != 0:
            sys.exit(f"[FATAL] exporter failed (exit {rc}). "
                     "Check you are in a PyTorch env and CLAIR3_HOME/--clair3 is set.")

        if registry_path.is_file():
            registry = yaml.safe_load(registry_path.read_text()) or {}
        else:
            registry = {}
        registry.setdefault("server_version", "0.0.1")
        registry.setdefault("models", [])

        want = STAGES if args.stage == "both" else (args.stage,)
        results = []
        for stage in want:
            src_ts = Path(tmp) / f"{stage}_ont_hac_v500.ts"
            src_meta = Path(tmp) / f"{stage}_ont_hac_v500.ts.meta.json"
            if not src_ts.is_file():
                print(f"[convert] note: {stage} not produced (no {stage}.pt in --model-dir?) — skipping")
                continue
            model_id = f"{stage}_{args.name}"
            dst_ts = models_dir / f"{model_id}.ts"
            shutil.copyfile(src_ts, dst_ts)
            ts_sha = sha256_file(dst_ts)
            meta = json.loads(src_meta.read_text())
            ts_rel = os.path.relpath(dst_ts, registry_dir)
            entry = build_entry(meta, model_id, args.platform, stage, args.accuracy, ts_rel, ts_sha)
            action = upsert(registry, entry)
            results.append((model_id, action, dst_ts, ts_sha))

        if not results:
            sys.exit("[FATAL] nothing exported — check --model-dir contains pileup.pt / full_alignment.pt")

        registry_path.write_text(yaml.safe_dump(registry, sort_keys=False, default_flow_style=False))

    print("\n[convert] done:")
    for model_id, action, dst_ts, ts_sha in results:
        print(f"  {action:7} {model_id}  ->  {dst_ts}  (sha256 {ts_sha[:16]}…)")
    print(f"  registry: {registry_path}")
    print("\nStart the server:")
    print(f"  clair3-connect-server --model-registry {registry_path} --device cpu --non-interactive")
    print("Or manage models interactively:  clair3-connect-server --interactive")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
