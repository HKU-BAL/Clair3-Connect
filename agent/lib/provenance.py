#!/usr/bin/env python3
"""Provenance manifest writer — stamps every results file with the context
needed to trust + reproduce it.

Usage (from any harness script):
    python lib/provenance.py \
        --out results/e1_parity/manifest.json \
        --experiment e1_parity \
        --extra key=value key2=value2 \
        --input-sha BAM=/path/to.bam

Captures: hostname, OS, CPU count, GPUs (nvidia-smi), git sha + dirty flag,
UTC timestamp, python/torch versions if importable, plus any caller-supplied
key=value pairs and sha256 of named input files.

Results files themselves are written by each harness; this only writes the
sidecar manifest so figures/ never has to guess where data came from.
"""
import argparse
import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
from datetime import datetime, timezone

def _run(cmd):
    try:
        return subprocess.check_output(cmd, shell=True, text=True,
                                       stderr=subprocess.DEVNULL).strip()
    except Exception:
        return ""

def _gpus():
    if not shutil.which("nvidia-smi"):
        return []
    out = _run("nvidia-smi --query-gpu=index,name,memory.total,driver_version "
               "--format=csv,noheader")
    return [line.strip() for line in out.splitlines() if line.strip()]

def _sha256(path, limit_mb=None):
    """sha256 of a file. For huge BAMs, hashing the whole file is slow, so
    callers can pass limit_mb to hash only a prefix (recorded as such)."""
    h = hashlib.sha256()
    n = 0
    cap = None if limit_mb is None else limit_mb * 1024 * 1024
    try:
        with open(path, "rb") as f:
            while True:
                chunk = f.read(1 << 20)
                if not chunk:
                    break
                if cap is not None and n + len(chunk) > cap:
                    h.update(chunk[: cap - n])
                    return h.hexdigest() + f"  (prefix:{limit_mb}MB)"
                h.update(chunk)
                n += len(chunk)
        return h.hexdigest()
    except Exception as e:
        return f"<unreadable: {e}>"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--experiment", required=True)
    ap.add_argument("--extra", nargs="*", default=[],
                    help="key=value pairs to record verbatim")
    ap.add_argument("--input-sha", nargs="*", default=[],
                    help="NAME=path pairs to sha256 (full file)")
    ap.add_argument("--input-sha-prefix", nargs="*", default=[],
                    help="NAME=path pairs to sha256 (first 64MB only, for huge BAMs)")
    ap.add_argument("--git-root", default=os.environ.get("CLAIR3_CONNECT_ROOT", "."))
    args = ap.parse_args()

    git_sha = _run(f"cd {args.git_root} && git rev-parse HEAD")
    git_dirty = bool(_run(f"cd {args.git_root} && git status --porcelain"))

    versions = {}
    try:
        import torch
        versions["torch"] = torch.__version__
        versions["torch_cuda"] = torch.version.cuda
    except Exception:
        pass
    versions["python"] = sys.version.split()[0]

    manifest = {
        "experiment": args.experiment,
        "utc": datetime.now(timezone.utc).isoformat(),
        "host": platform.node(),
        "os": platform.platform(),
        "cpu_count": os.cpu_count(),
        "gpus": _gpus(),
        "git_sha": git_sha,
        "git_dirty": git_dirty,
        "versions": versions,
        "extra": dict(kv.split("=", 1) for kv in args.extra if "=" in kv),
        "input_sha256": {},
    }
    for kv in args.input_sha:
        if "=" in kv:
            name, path = kv.split("=", 1)
            manifest["input_sha256"][name] = {"path": path, "sha256": _sha256(path)}
    for kv in args.input_sha_prefix:
        if "=" in kv:
            name, path = kv.split("=", 1)
            manifest["input_sha256"][name] = {"path": path,
                                              "sha256": _sha256(path, limit_mb=64)}

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"[provenance] wrote {args.out}")

if __name__ == "__main__":
    main()
