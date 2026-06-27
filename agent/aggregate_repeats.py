#!/usr/bin/env python3

import glob
import json
import os
import statistics as st
import sys

ARMS = ["tools", "skill", "baseline"]
METRICS = ["turns", "tool_calls", "tokens", "wall_s"]

def msd(xs):
    """(mean, sample-sd, n). sd is None for n<2."""
    n = len(xs)
    if n == 0:
        return (None, None, 0)
    m = st.mean(xs)
    s = st.stdev(xs) if n > 1 else None
    return (round(m, 2), (round(s, 2) if s is not None else None), n)

def main():
    base = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
        os.environ.get("RESULTS", "."), "e4_agent", "repeats")
    repdirs = sorted(d for d in glob.glob(os.path.join(base, "rep*")) if os.path.isdir(d))
    if not repdirs:
        sys.exit(f"no rep* dirs under {base}")

    rows = {a: [] for a in ARMS}
    n_reps_seen = {a: 0 for a in ARMS}
    for rd in repdirs:
        for arm in ARMS:
            f = os.path.join(rd, f"{arm}.only_transcripts.json")
            if not os.path.isfile(f):
                continue
            n_reps_seen[arm] += 1
            for r in json.load(open(f)):
                rows[arm].append({
                    "rep": os.path.basename(rd),
                    "task_id": r["task_id"],
                    "success": bool(r["success"]),
                    "turns": r["turns"],
                    "tool_calls": r["tool_calls"],
                    "tokens": (r.get("prompt_tokens", 0) + r.get("completion_tokens", 0)),
                    "prompt_tokens": r.get("prompt_tokens", 0),
                    "completion_tokens": r.get("completion_tokens", 0),
                    "wall_s": r["wall_s"],
                    "forced_answer": bool(r.get("forced_answer")),
                })

    summary = {"base": base, "n_repeats_per_arm": n_reps_seen, "arms": {}}
    for arm in ARMS:
        rs = rows[arm]
        if not rs:
            continue
        tasks = sorted({r["task_id"] for r in rs})
        per_sample = {}
        for tid in tasks:
            sub = [r for r in rs if r["task_id"] == tid]
            per_sample[tid] = {
                "n": len(sub),
                "success_rate": round(sum(r["success"] for r in sub) / len(sub), 3),
                **{mk: msd([r[mk] for r in sub]) for mk in METRICS},
            }
        pooled = {
            "n_runs": len(rs),
            "success_rate": round(sum(r["success"] for r in rs) / len(rs), 3),
            "forced_answers": sum(r["forced_answer"] for r in rs),
            **{mk: msd([r[mk] for r in rs]) for mk in METRICS},
        }

        toks = [r["tokens"] for r in rs]
        cv = round(st.stdev(toks) / st.mean(toks), 3) if len(toks) > 1 and st.mean(toks) else None
        pooled["tokens_cv"] = cv
        summary["arms"][arm] = {"pooled": pooled, "per_sample": per_sample}

    if "tools" in summary["arms"]:
        t_tok = summary["arms"]["tools"]["pooled"]["tokens"][0]
        ratios = {}
        for arm in ARMS:
            if arm in summary["arms"] and t_tok:
                ratios[arm] = round(summary["arms"][arm]["pooled"]["tokens"][0] / t_tok, 2)
        summary["token_ratio_vs_tools"] = ratios

    out = os.path.join(base, "summary.json")
    json.dump(summary, open(out, "w"), indent=2)

    print(f"\nE4 APOE repeat aggregate  ({base})")
    print(f"repeats per arm: {n_reps_seen}\n")
    hdr = f"{'arm':9} {'n':>3} {'succ':>5} {'turns':>13} {'tool_calls':>13} {'tokens':>17} {'wall_s':>13} {'tokCV':>6} {'×tools':>7}"
    print(hdr)
    print("-" * len(hdr))
    for arm in ARMS:
        if arm not in summary["arms"]:
            continue
        p = summary["arms"][arm]["pooled"]
        def fmt(mk):
            m, s, _ = p[mk]
            return f"{m}±{s}" if s is not None else f"{m}"
        ratio = summary.get("token_ratio_vs_tools", {}).get(arm, "")
        print(f"{arm:9} {p['n_runs']:>3} {p['success_rate']:>5} "
              f"{fmt('turns'):>13} {fmt('tool_calls'):>13} {fmt('tokens'):>17} "
              f"{fmt('wall_s'):>13} {str(p['tokens_cv']):>6} {str(ratio):>7}")
    print("\nper-sample token mean±sd:")
    for arm in ARMS:
        if arm not in summary["arms"]:
            continue
        for tid, d in summary["arms"][arm]["per_sample"].items():
            m, s, n = d["tokens"]
            print(f"  {arm:9} {tid:28} n={n} succ={d['success_rate']} tokens={m}±{s}")
    print(f"\nwrote {out}")

if __name__ == "__main__":
    main()
