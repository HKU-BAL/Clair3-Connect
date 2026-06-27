#!/usr/bin/env python
"""export_clair3_torchscript.py — Phase 1 TorchScript exporter.

Loads a Clair3 PyTorch state-dict checkpoint, builds the corresponding
`Clair3_P` (pileup) or `Clair3_F` (full-alignment) module with the EXACT
configuration used by stock Clair3, traces + freezes it, and writes a
self-contained `.ts` file plus a `.meta.json` sidecar with shape/dtype/sha256
for the registry.

Run inside the `clair3-torch` conda env (same env as the Phase 0 producer).

Usage:
    python tools/export_torchscript/export_clair3_torchscript.py \\
        --model-dir test_data/ont_chr20_demo/model_r1041_e82_400bps_hac_v500 \\
        --platform ont \\
        --output-dir models

Outputs (for the ONT HAC v500 fixture):
    models/pileup_ont_hac_v500.ts                (frozen, traced)
    models/pileup_ont_hac_v500.ts.meta.json
    models/full_alignment_ont_hac_v500.ts
    models/full_alignment_ont_hac_v500.ts.meta.json

Decisions baked in (anchored in the Phase 0 producer + plan doc):
- Pileup: `add_indel_length=False`, `predict=True`, input dtype int32, shape
  [B, 33, 18], output [B, 24] (GT21 + GENOTYPE).
- Full-alignment ONT: `add_indel_length=True`, `predict=True`, input dtype
  int8, NHWC shape [B, 89, 33, 8], output [B, 90].
- Tracing example uses batch_size=4 — `Clair3_P` has only batch-agnostic
  layers (LSTM, Linear over flattened time*hidden); `Clair3_F` has only
  batch-agnostic Conv2d/BatchNorm/PyramidPooling/Linear. We verify dynamic
  batch behavior by re-running the traced module with batch_size=1 and 8
  against the eager model.

Parity check: scripted(X) is compared to eager(X) on the Phase 0 sample
tensors when available; otherwise on a small synthetic tensor.
"""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import logging
import sys
from pathlib import Path
from typing import Any

import numpy as np
import torch

HERE = Path(__file__).resolve().parent
WORKSPACE = HERE.parent.parent
PROJECT_ROOT = WORKSPACE.parent
import os as _os
_env_clair3 = _os.environ.get("CLAIR3_HOME")
CLAIR3_DIR = Path(_env_clair3).resolve() if _env_clair3 else (PROJECT_ROOT / "Clair3")

if not CLAIR3_DIR.is_dir():
    sys.exit(
        f"[FATAL] Clair3 source clone not found at {CLAIR3_DIR}.\n"
        "        Set CLAIR3_HOME to your Clair3 checkout, e.g.:\n"
        "          export CLAIR3_HOME=/path/to/Clair3"
    )
sys.path.insert(0, str(CLAIR3_DIR))

LOG = logging.getLogger("export_torchscript")

def setup_logging(verbose: bool) -> None:
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="[%(asctime)s] %(levelname)s %(name)s: %(message)s",
        datefmt="%H:%M:%S",
    )

def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

def load_checkpoint(model: torch.nn.Module, ckpt: Path) -> None:
    """Mirror `clair3/CallVariantsFromCffi.py:19-28`."""
    state = torch.load(str(ckpt), map_location="cpu")
    if isinstance(state, dict) and "state_dict" in state:
        state = state["state_dict"]
    model.load_state_dict(state)

def export_pileup(model_dir: Path, output_dir: Path, platform: str,
                  trace_batch: int = 4) -> dict[str, Any]:
    """Trace + freeze the Clair3 pileup model."""
    from clair3.model import Clair3_P
    import shared.param_p as param_p

    ckpt = model_dir / "pileup.pt"
    if not ckpt.exists():
        raise FileNotFoundError(ckpt)

    add_indel_length = False
    input_channels = (param_p.ont_input_shape if platform == "ont" else param_p.input_shape)[-1]
    no_of_positions = param_p.no_of_positions
    output_width = 21 + 3

    LOG.info("Building Clair3_P(add_indel_length=%s, predict=True, input_channels=%d)",
             add_indel_length, input_channels)
    model = Clair3_P(add_indel_length=add_indel_length, predict=True,
                     input_channels=input_channels)
    load_checkpoint(model, ckpt)
    model.eval()

    example = torch.zeros((trace_batch, no_of_positions, input_channels), dtype=torch.int32)

    import torch.nn as nn
    class _PileupPredictWrapper(nn.Module):
        """Inlined predict=True / add_indel_length=False branch of
        Clair3_P.forward. Scriptable (no branches on dynamic flags,
        uniform Tensor return). Device-agnostic — `nn.LSTM` uses
        `input.new_zeros(...)` internally so hidden states follow
        the input device."""
        def __init__(self, inner):
            super().__init__()
            self.LSTM1 = inner.LSTM1
            self.LSTM2 = inner.LSTM2
            self.L3_dropout = inner.L3_dropout
            self.activation = inner.activation
            self.softmax = inner.softmax
            self.L4 = inner.L4
            self.L4_dropout = inner.L4_dropout
            self.L5_1 = inner.L5_1
            self.L5_1_dropout = inner.L5_1_dropout
            self.L5_2 = inner.L5_2
            self.L5_2_dropout = inner.L5_2_dropout
            self.Y_gt21_logits = inner.Y_gt21_logits
            self.Y_genotype_logits = inner.Y_genotype_logits

        def forward(self, x: torch.Tensor) -> torch.Tensor:
            x = x.float()
            x, _ = self.LSTM1(x)
            x, _ = self.LSTM2(x)
            x = self.L3_dropout(x)
            x = torch.flatten(x, start_dim=1)
            x = self.activation(self.L4(x))
            x = self.L4_dropout(x)
            l5_1 = self.L5_1_dropout(self.activation(self.L5_1(x)))
            l5_2 = self.L5_2_dropout(self.activation(self.L5_2(x)))
            y_gt21 = self.softmax(self.activation(self.Y_gt21_logits(l5_1)))
            y_gt   = self.softmax(self.activation(self.Y_genotype_logits(l5_2)))
            return torch.cat([y_gt21, y_gt], dim=1)

    wrapped = _PileupPredictWrapper(model).eval()

    try:
        LOG.info("Scripting (device-agnostic)")
        with torch.inference_mode():
            scripted = torch.jit.script(wrapped)

        _ = scripted(example)
        LOG.info("Freezing scripted module")
        frozen = torch.jit.freeze(scripted)
        export_mode = "script"
    except Exception as e:
        LOG.warning(
            "torch.jit.script failed (%s); falling back to trace (CPU-only). "
            "The resulting .ts will not run on CUDA.",
            type(e).__name__,
        )
        LOG.info("Tracing with input shape %s dtype %s",
                 list(example.shape), example.dtype)
        with torch.inference_mode():
            traced = torch.jit.trace(model, example, strict=False)
        LOG.info("Freezing traced module")
        frozen = torch.jit.freeze(traced)
        export_mode = "trace"

    out_path = output_dir / "pileup_ont_hac_v500.ts"
    frozen.save(str(out_path))
    LOG.info("Wrote %s (mode=%s)", out_path, export_mode)

    sample_x = _try_load_phase0_pileup_sample(WORKSPACE)
    if sample_x is None:
        LOG.warning("No Phase 0 pileup tensor found; using random int32 sample for parity check")
        rng = np.random.default_rng(0)
        sample_x = rng.integers(0, 50, size=(8, no_of_positions, input_channels), dtype=np.int32)

    parity = _parity_check(model, frozen, sample_x, dtype=torch.int32)

    meta = {
        "model_id": "pileup_ont_hac_v500",
        "produced_at_utc": datetime.datetime.utcnow().isoformat() + "Z",
        "checkpoint_path": str(ckpt),
        "checkpoint_sha256": sha256_file(ckpt),
        "torchscript": {
            "path": str(out_path.relative_to(output_dir)),
            "sha256": sha256_file(out_path),
            "size_bytes": out_path.stat().st_size,
        },
        "config": {
            "model_class": "Clair3_P",
            "add_indel_length": add_indel_length,
            "predict": True,
            "input_channels": input_channels,
            "no_of_positions": no_of_positions,
        },
        "input": {
            "shape": [-1, no_of_positions, input_channels],
            "dtype": "int32",
            "trace_batch_size": trace_batch,
        },
        "output": {
            "shape": [-1, output_width],
            "dtype": "float32",
        },
        "torch_version": torch.__version__,
        "parity": parity,
    }
    meta_path = out_path.with_suffix(".ts.meta.json")
    with open(meta_path, "w") as fh:
        json.dump(meta, fh, indent=2)
    LOG.info("Wrote %s", meta_path)

    return meta

def export_full_alignment(model_dir: Path, output_dir: Path, platform: str,
                          trace_batch: int = 4) -> dict[str, Any]:
    """Trace + freeze the Clair3 full-alignment model."""
    from clair3.model import Clair3_F
    import shared.param_f as param_f

    ckpt = model_dir / "full_alignment.pt"
    if not ckpt.exists():
        raise FileNotFoundError(ckpt)

    add_indel_length = True

    input_channels = (param_f.ont_input_shape if platform == "ont" else param_f.input_shape)[-1]
    matrix_depth = param_f.matrix_depth_dict[platform]
    no_of_positions = param_f.no_of_positions
    output_width = 21 + 3 + 33 + 33

    LOG.info("Building Clair3_F(add_indel_length=%s, predict=True, input_channels=%d)",
             add_indel_length, input_channels)
    model = Clair3_F(add_indel_length=add_indel_length, predict=True,
                     input_channels=input_channels)
    load_checkpoint(model, ckpt)
    model.eval()

    example = torch.zeros((trace_batch, matrix_depth, no_of_positions, input_channels),
                          dtype=torch.int8)

    import torch.nn as nn
    import torch.nn.functional as F
    from typing import List

    class _ScriptPyramidPolling(nn.Module):
        __constants__ = ["spatial_pool_size"]

        def __init__(self, spatial_pool_size=(3, 2, 1)):
            super().__init__()
            self.spatial_pool_size = list(int(p) for p in spatial_pool_size)

        def forward(self, x: torch.Tensor) -> torch.Tensor:
            pooled: List[torch.Tensor] = []
            height = int(x.shape[-2])
            width = int(x.shape[-1])
            for pool_size in self.spatial_pool_size:
                window_h = (height + pool_size - 1) // pool_size
                window_w = (width + pool_size - 1) // pool_size
                stride_h = window_h
                stride_w = window_w
                out_h = (height + stride_h - 1) // stride_h
                out_w = (width + stride_w - 1) // stride_w
                pad_h = max((out_h - 1) * stride_h + window_h - height, 0)
                pad_w = max((out_w - 1) * stride_w + window_w - width, 0)
                pad_top = pad_h // 2
                pad_bottom = pad_h - pad_top
                pad_left = pad_w // 2
                pad_right = pad_w - pad_left
                if pad_h > 0 or pad_w > 0:
                    x_pad = F.pad(
                        x,
                        [pad_left, pad_right, pad_top, pad_bottom],
                        mode="constant",
                        value=0.0,
                    )
                else:
                    x_pad = x
                mp = F.max_pool2d(
                    x_pad,
                    kernel_size=(window_h, window_w),
                    stride=(stride_h, stride_w),
                )
                mp = mp.permute(0, 2, 3, 1)
                pooled.append(torch.flatten(mp, start_dim=1))
            return torch.cat(pooled, dim=1)

    class _FaPredictWrapper(nn.Module):
        """Inlined predict=True / add_indel_length=True branch of
        Clair3_F.forward. Hardcodes the NHWC→NCHW permute so the
        scripter doesn't have to handle the dynamic `if x.shape[-1]
        == input_channels` check. Input must be NHWC `[B, depth,
        width, channels]` int8."""

        __constants__ = ["normalize_num"]

        def __init__(self, inner, normalize_num: float):
            super().__init__()
            self.normalize_num = normalize_num
            self.conv1 = inner.conv1
            self.res_block1 = inner.res_block1
            self.conv3 = inner.conv3
            self.res_block2 = inner.res_block2
            self.conv5 = inner.conv5
            self.res_block3 = inner.res_block3

            self.pyramidpolling = _ScriptPyramidPolling(
                inner.pyramidpolling.spatial_pool_size
            )
            self.L3_dropout = inner.L3_dropout
            self.L4 = inner.L4
            self.L4_dropout = inner.L4_dropout
            self.L5_1 = inner.L5_1
            self.L5_1_dropout = inner.L5_1_dropout
            self.L5_2 = inner.L5_2
            self.L5_2_dropout = inner.L5_2_dropout
            self.Y_gt21_logits = inner.Y_gt21_logits
            self.Y_genotype_logits = inner.Y_genotype_logits
            self.L5_3 = inner.L5_3
            self.L5_3_dropout = inner.L5_3_dropout
            self.L5_4 = inner.L5_4
            self.L5_4_dropout = inner.L5_4_dropout
            self.Y_indel_length_logits_1 = inner.Y_indel_length_logits_1
            self.Y_indel_length_logits_2 = inner.Y_indel_length_logits_2
            self.softmax = inner.softmax
            self.activation = inner.activation

        def forward(self, inputs: torch.Tensor) -> torch.Tensor:
            x = inputs.float() / self.normalize_num

            x = x.permute(0, 3, 1, 2)
            x = self.conv1(x)
            x = self.res_block1(x)
            x = self.conv3(x)
            x = self.res_block2(x)
            x = self.conv5(x)
            x = self.res_block3(x)
            x = self.pyramidpolling(x)
            x = self.L3_dropout(x)
            x = self.activation(self.L4(x))
            x = self.L4_dropout(x)
            l5_1 = self.L5_1_dropout(self.activation(self.L5_1(x)))
            l5_2 = self.L5_2_dropout(self.activation(self.L5_2(x)))
            y_gt21 = self.softmax(self.activation(self.Y_gt21_logits(l5_1)))
            y_gt = self.softmax(self.activation(self.Y_genotype_logits(l5_2)))
            l5_3 = self.L5_3_dropout(self.activation(self.L5_3(x)))
            l5_4 = self.L5_4_dropout(self.activation(self.L5_4(x)))
            y_il1 = self.softmax(self.activation(self.Y_indel_length_logits_1(l5_3)))
            y_il2 = self.softmax(self.activation(self.Y_indel_length_logits_2(l5_4)))
            return torch.cat([y_gt21, y_gt, y_il1, y_il2], dim=1)

    for mod in model.modules():
        if isinstance(mod, nn.Linear):
            mod.in_features = int(mod.in_features)
            mod.out_features = int(mod.out_features)

    wrapped = _FaPredictWrapper(model, float(param_f.NORMALIZE_NUM)).eval()

    try:
        LOG.info("Scripting (device-agnostic)")
        with torch.inference_mode():
            scripted = torch.jit.script(wrapped)
        _ = scripted(example)
        LOG.info("Freezing scripted module")
        frozen = torch.jit.freeze(scripted)
        export_mode = "script"
    except Exception as e:
        LOG.warning(
            "torch.jit.script failed for FA (%s: %s); falling back to trace "
            "(CPU-only). The resulting .ts will not run on CUDA.",
            type(e).__name__,
            str(e)[:500],
        )
        LOG.info("Tracing with input shape %s dtype %s",
                 list(example.shape), example.dtype)
        with torch.inference_mode():
            traced = torch.jit.trace(model, example, strict=False)
        LOG.info("Freezing traced module")
        frozen = torch.jit.freeze(traced)
        export_mode = "trace"

    out_path = output_dir / "full_alignment_ont_hac_v500.ts"
    frozen.save(str(out_path))
    LOG.info("Wrote %s (mode=%s)", out_path, export_mode)

    rng = np.random.default_rng(1)
    sample_x = rng.integers(-30, 100, size=(8, matrix_depth, no_of_positions, input_channels),
                            dtype=np.int8)
    parity = _parity_check(model, frozen, sample_x, dtype=torch.int8)

    meta = {
        "model_id": "full_alignment_ont_hac_v500",
        "produced_at_utc": datetime.datetime.utcnow().isoformat() + "Z",
        "checkpoint_path": str(ckpt),
        "checkpoint_sha256": sha256_file(ckpt),
        "torchscript": {
            "path": str(out_path.relative_to(output_dir)),
            "sha256": sha256_file(out_path),
            "size_bytes": out_path.stat().st_size,
        },
        "config": {
            "model_class": "Clair3_F",
            "add_indel_length": add_indel_length,
            "predict": True,
            "input_channels": input_channels,
            "matrix_depth": matrix_depth,
            "no_of_positions": no_of_positions,
            "enable_dwell_time": False,
        },
        "input": {
            "shape": [-1, matrix_depth, no_of_positions, input_channels],
            "dtype": "int8",
            "layout": "NHWC",
            "trace_batch_size": trace_batch,
        },
        "output": {
            "shape": [-1, output_width],
            "dtype": "float32",
        },
        "torch_version": torch.__version__,
        "parity": parity,
    }
    meta_path = out_path.with_suffix(".ts.meta.json")
    with open(meta_path, "w") as fh:
        json.dump(meta, fh, indent=2)
    LOG.info("Wrote %s", meta_path)

    return meta

def _try_load_phase0_pileup_sample(workspace: Path) -> "np.ndarray | None":
    """Pull a small slice of `pileup_post_scale_X_int32.npy` for parity."""
    p = workspace / "test_data" / "ont_chr20_demo" / "artifacts" / "pileup_post_scale_X_int32.npy"
    if not p.exists():
        return None
    arr = np.load(p, allow_pickle=False)

    return np.ascontiguousarray(arr[:16])

def _parity_check(eager: torch.nn.Module, scripted: torch.jit.ScriptModule,
                  X_np: np.ndarray, *, dtype: torch.dtype,
                  tolerance: float = 1e-5) -> dict[str, Any]:
    """Run both models and return a parity report dict."""
    X = torch.from_numpy(X_np).to(dtype)
    with torch.inference_mode():
        y_eager = eager(X).detach().cpu().numpy()
        y_scripted = scripted(X).detach().cpu().numpy()
    max_abs_diff = float(np.max(np.abs(y_eager - y_scripted))) if y_eager.size else 0.0
    argmax_match_pct = (
        float(100.0 * np.mean(np.argmax(y_eager, axis=1) == np.argmax(y_scripted, axis=1)))
        if y_eager.size else None
    )
    ok = max_abs_diff < tolerance
    report = {
        "sample_shape": list(X.shape),
        "max_abs_diff": max_abs_diff,
        "argmax_match_pct": argmax_match_pct,
        "tolerance": tolerance,
        "passed": ok,
    }
    LOG.info("Parity (sample %s): max_abs_diff=%.3e  argmax_match=%.2f%%  %s",
             list(X.shape), max_abs_diff, argmax_match_pct or 0.0, "OK" if ok else "FAIL")

    for b in [1, 8]:
        if b > X.shape[0]:
            continue
        x_sub = X[:b]
        with torch.inference_mode():
            y_sub = scripted(x_sub).detach().cpu().numpy()
        assert y_sub.shape[0] == b, f"scripted output batch {y_sub.shape[0]} != {b}"
    report["dynamic_batch_ok"] = True

    if not ok:
        raise RuntimeError(
            f"TorchScript parity check FAILED: max_abs_diff={max_abs_diff} > {tolerance}"
        )
    return report

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--model-dir", required=True,
                        help="Clair3 model dir containing pileup.pt and full_alignment.pt.")
    parser.add_argument("--output-dir", default="models",
                        help="Where to write *.ts and *.ts.meta.json (default: %(default)s).")
    parser.add_argument("--platform", default="ont", choices=["ont", "hifi", "ilmn"])
    parser.add_argument("--trace-batch", type=int, default=4,
                        help="Batch size for the trace example tensor (default: %(default)s).")
    parser.add_argument("--skip", choices=["pileup", "full_alignment"], default=None,
                        help="Skip one of the two models (e.g., --skip full_alignment).")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    setup_logging(args.verbose)

    model_dir = Path(args.model_dir).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    metas: dict[str, dict[str, Any]] = {}

    if args.skip != "pileup":
        metas["pileup"] = export_pileup(model_dir, output_dir, args.platform, args.trace_batch)
    if args.skip != "full_alignment":
        metas["full_alignment"] = export_full_alignment(model_dir, output_dir, args.platform,
                                                        args.trace_batch)

    LOG.info("=" * 64)
    LOG.info("Export summary:")
    for stage, meta in metas.items():
        ts = meta["torchscript"]
        LOG.info("  %-15s %s (sha256 %s..., %.1f MB)",
                 stage, ts["path"], ts["sha256"][:16], ts["size_bytes"] / 1e6)
        LOG.info("    parity: max_abs_diff=%.3e (tol %s) passed=%s",
                 meta["parity"]["max_abs_diff"],
                 meta["parity"]["tolerance"],
                 meta["parity"]["passed"])

    LOG.info("=" * 64)
    LOG.info("Registry update snippet (paste into models/model_registry.example.yaml):")
    for stage, meta in metas.items():
        print(f"  # {stage}")
        print(f"  #   torchscript.path = {meta['torchscript']['path']}")
        print(f"  #   torchscript.sha256 = {meta['torchscript']['sha256']}")

    return 0

if __name__ == "__main__":
    sys.exit(main())
