#!/usr/bin/env python3
import argparse
import json
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import clair3_connect_tools as T

def run_shell(command: str, timeout: int = 1800) -> dict:
    try:
        p = subprocess.run(command, shell=True, capture_output=True, text=True,
                           timeout=timeout, executable="/bin/bash")
        return {"returncode": p.returncode, "stdout": p.stdout[-4000:], "stderr": p.stderr[-1500:]}
    except subprocess.TimeoutExpired:
        return {"error": "timeout", "timeout_s": timeout}

BASELINE_TOOLS = [{"type": "function", "function": {
    "name": "run_shell",
    "description": "Run one bash command on the analysis host; returns returncode/stdout/stderr.",
    "parameters": {"type": "object", "properties": {
        "command": {"type": "string"}}, "required": ["command"]}}}]
BASELINE_DISPATCH = {"run_shell": lambda command: run_shell(command)}

SYS_TOOLS = (
    "You are a genomics analysis assistant. You can call germline variants and "
    "inspect results ONLY through the provided Clair3-connect tools. Datasets are "
    "referenced by (sample, platform); you never handle file paths. Think step "
    "by step, call tools as needed, and finish with a single final line of the "
    "exact form 'ANSWER: <value>' (a bare number when a count/fraction is asked)."
)

def _apoe_env_lines(clinical=False):
    """List every registered APOE clinical sample (from apoe_truth.json) with its
    BAM path, so the shell-based (baseline/skill) agents can LOCATE the data
    without scanning the filesystem. Platform is intentionally NOT given — the
    agent must read the BAM @PG/@RG header and pick the matching model (this is
    the model-selection part of the test). `clinical` adjusts the heading: on an
    APOE-only run these are the ONLY samples, not "additional" ones."""
    import json as _json
    try:
        d = _json.load(open(os.environ.get("APOE_TRUTH_JSON", "")))
    except Exception:
        return ""
    head = "Registered clinical samples" if clinical else "Additional registered clinical samples"
    out = [head + " (human WGS; determine each one's platform yourself from its "
           "BAM @PG/@RG header, then pick the matching model/ref):"]
    for s, info in sorted(d.items()):
        out.append(f"  {s} bam={info.get('bam','')} ref={os.environ.get('REF_ILMN','')}")
    return "\n".join(out) + "\n"

def _local_env_block(clinical=False):
    """Local execution env shared by the baseline + skill conditions (this host
    runs stock Clair3 via conda, not the Docker image the skill assumes).

    The block injects ONLY the clinically-relevant context: the run_clair3
    template, the model dirs (so the agent must still SELECT the matching model
    from its BAM header), and the registered APOE samples with the single GRCh38
    reference they use."""
    c = os.environ.get("CLAIR3_ROOT", "")
    head = (
        "Local execution environment (NO Docker image here; run stock Clair3 "
        "locally):\n"
        f"  conda run -n clair3-torch python {c}/run_clair3.py --bam_fn BAM "
        "--ref_fn REF --model_path MODELDIR --platform {ont|hifi|ilmn} --threads 8 "
        f"--output OUTDIR --bed_fn REGION.bed --device cpu --longphase {c}/longphase "
        "--sample_name SAMPLE\n"
        "Model dirs (each holds pileup.pt + full_alignment.pt):\n"
        f"  ont={os.environ.get('CKPT_ONT','')}\n  hifi={os.environ.get('CKPT_HIFI','')}\n"
        f"  ilmn={os.environ.get('CKPT_ILMN','')}\n"
    )
    return (
        head
        + _apoe_env_lines(clinical=True)
        + "Output is OUTDIR/merge_output.vcf.gz (bgzipped). Build any region BED "
          "yourself.\n"
        f"Read genotypes with bcftools ({os.environ.get('BCFTOOLS','bcftools')})."
    )

def sys_baseline(clinical=False):
    return (
        "You are a genomics analysis assistant on a Linux host with a single tool: "
        "run_shell, and a pointer to the stock (monolithic) Clair3 CLI.\n"
        + _local_env_block(clinical)
        + "\nFinish with a single final line 'ANSWER: <value>'."
    )

def sys_skill(clinical=False):
    """Skill condition: the agent is given the Claude-Code-style Clair skill
    (HKU-BAL/Clair-skills SKILL.md) + the same local exec env, and may read the
    skill's references/ and templates/ via run_shell. Mirrors how an LLM agent
    consumes a packaged skill (large instruction doc + bundled scripts + shell)."""
    skill_dir = os.environ.get("CLAIR_SKILL_DIR", "")
    try:
        skill_md = open(os.path.join(skill_dir, "SKILL.md")).read()
    except Exception:
        skill_md = "(SKILL.md unavailable)"
    return (
        "You are a genomics analysis assistant equipped with the SKILL below "
        "(read and follow it) and a single tool: run_shell.\n\n"
        "===== BEGIN SKILL.md =====\n" + skill_md + "\n===== END SKILL.md =====\n\n"
        f"The skill's reference docs live in {skill_dir}/references/*.md and command "
        f"templates in {skill_dir}/templates/*.sh — read any with run_shell (`cat`). "
        + _local_env_block(clinical)
        + "\nFinish with a single final line 'ANSWER: <value>'."
    )

def _assistant_to_dict(msg):
    d = {"role": "assistant", "content": msg.content or ""}

    rc = getattr(msg, "reasoning_content", None)
    if rc is None and getattr(msg, "model_extra", None):
        rc = msg.model_extra.get("reasoning_content")
    if rc:
        d["reasoning_content"] = rc
    if msg.tool_calls:
        d["tool_calls"] = [{
            "id": tc.id, "type": "function",
            "function": {"name": tc.function.name, "arguments": tc.function.arguments},
        } for tc in msg.tool_calls]
    return d

def _force_final_answer(client, model, messages, reason):
    """The agent hit its step/time budget while still mid-analysis (its last
    action was a tool call, so no free-text 'ANSWER:' was ever produced). Ask
    it ONCE more, with tools disabled, to commit the best conclusion from the
    evidence ALREADY gathered.

    Why this is a fairness fix, not a thumb on the scale:
    the raw-CLI baseline needs many more shell turns than the structured-tools
    arm to run the same analysis, so a fixed turn cap can truncate a baseline
    run that has, in fact, already computed the answer (e.g. concordance=0.94,
    precision=0.997) — scoring a COMPLETED analysis as a null answer. Forcing a
    final answer measures "did the analysis reach the right conclusion" rather
    than "did the agent leave itself a spare turn to print it". It is applied
    IDENTICALLY to every arm and never injects any task-specific information; it
    is a no-op for the tools/skill arms, which answer well within budget. The
    cost of needing the nudge is still captured honestly in turns/tokens.
    Returns (answer_text, usage)."""
    msgs = messages + [{"role": "user", "content": (
        f"You have reached your analysis budget ({reason}). Do NOT call any more "
        "tools. Based ONLY on the evidence you have already gathered above, state "
        "your single best final conclusion now as one line: 'ANSWER: <value>'.")}]
    try:
        resp = client.chat.completions.create(
            model=model, messages=msgs, tools=tools_param(messages),
            tool_choice="none", temperature=0, timeout=120)
        return (resp.choices[0].message.content or ""), resp.usage
    except Exception as e:
        return f"(force-answer api-error: {e})", None

def tools_param(messages):
    """Some APIs reject tool_choice='none' unless a tools list is also present
    (the message history already contains tool_calls). Re-derive the right tool
    list from the system prompt rather than threading it through."""
    sysmsg = messages[0]["content"] if messages else ""
    return T.TOOLS if "Clair3-connect tools" in sysmsg else BASELINE_TOOLS

def run_task(client, model, task, mode, max_turns=None):
    if max_turns is None:

        max_turns = int(os.environ.get("E4_MAX_TURNS", "50"))

    clinical = task.get("match") == "diplotype"
    if mode == "tools":

        tools = T.TOOLS_CLINICAL if clinical else T.TOOLS
        dispatch, sys_prompt = T.DISPATCH, SYS_TOOLS
    elif mode == "skill":
        tools, dispatch, sys_prompt = BASELINE_TOOLS, BASELINE_DISPATCH, sys_skill(clinical)
    else:
        tools, dispatch, sys_prompt = BASELINE_TOOLS, BASELINE_DISPATCH, sys_baseline(clinical)
    messages = [{"role": "system", "content": sys_prompt},
                {"role": "user", "content": task["prompt"]}]
    m = {"turns": 0, "tool_calls": 0, "prompt_tokens": 0, "completion_tokens": 0,
         "tool_errors": 0, "tool_names": [], "forced_answer": False}
    final = ""
    answered = False
    budget_hit = None

    budget_s = float(os.environ.get("E4_TASK_BUDGET_S", "1800"))
    t0 = time.time()
    for _ in range(max_turns):
        if time.time() - t0 > budget_s:
            budget_hit = "wall-clock time budget"
            break
        m["turns"] += 1
        try:
            resp = client.chat.completions.create(
                model=model, messages=messages, tools=tools, temperature=0, timeout=120)
        except Exception as e:
            final = f"(api-error: {e})"
            answered = True
            break
        if resp.usage:
            m["prompt_tokens"] += resp.usage.prompt_tokens or 0
            m["completion_tokens"] += resp.usage.completion_tokens or 0
        msg = resp.choices[0].message
        messages.append(_assistant_to_dict(msg))
        if msg.tool_calls:
            for tc in msg.tool_calls:
                m["tool_calls"] += 1
                fn = tc.function.name
                m["tool_names"].append(fn)
                try:
                    args = json.loads(tc.function.arguments or "{}")
                    result = dispatch[fn](**args)
                except Exception as e:
                    result = {"error": str(e)}
                if isinstance(result, dict) and "error" in result:
                    m["tool_errors"] += 1
                messages.append({"role": "tool", "tool_call_id": tc.id,
                                 "content": json.dumps(result)[:6000]})
            continue
        final = msg.content or ""
        answered = True
        break
    if not answered and budget_hit is None:
        budget_hit = "step budget"

    if not answered and budget_hit:
        forced, usage = _force_final_answer(client, model, messages, budget_hit)
        if usage:
            m["prompt_tokens"] += usage.prompt_tokens or 0
            m["completion_tokens"] += usage.completion_tokens or 0
        final = forced
        m["forced_answer"] = True
        m["turns"] += 1
    m["wall_s"] = round(time.time() - t0, 2)
    ok = score(final, task)
    return {"task_id": task["id"], "mode": mode, "success": ok, "final": final[-400:],
            **m, "transcript": messages}

def score(final, task):
    mobj = re.search(r"ANSWER:\s*(.+)", final or "", re.IGNORECASE)
    ans = (mobj.group(1).strip() if mobj else (final or "").strip())
    exp = task["expected"]
    kind = task.get("match", "diplotype")
    try:
        if kind == "diplotype":

            def _alleles(s):
                s = s.lower().replace("ε", "e")
                m = re.search(r"(?:e)?([234])\s*/\s*(?:e)?([234])", s)
                return sorted([m.group(1), m.group(2)]) if m else []
            return _alleles(ans) == _alleles(str(exp)) and _alleles(ans) != []
    except Exception:
        return False
    return False

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tasks", default=os.path.join(HERE, "tasks.json"))
    ap.add_argument("--mode", choices=["tools", "skill", "baseline", "both", "all"], default="all")
    ap.add_argument("--out-dir", default=os.environ.get("RESULTS", ".") + "/e4_agent")
    ap.add_argument("--model", default=os.environ.get("DEEPSEEK_MODEL", "deepseek-chat"))
    ap.add_argument("--only", default=os.environ.get("E4_ONLY", ""),
                    help="comma-separated task ids to run (others skipped); empty = all")
    args = ap.parse_args()

    from openai import OpenAI
    key = os.environ.get("DEEPSEEK_API_KEY")
    if not key:
        sys.exit("DEEPSEEK_API_KEY not set — export it before running E4.")
    client = OpenAI(api_key=key, base_url=os.environ.get("DEEPSEEK_BASE_URL", "https://api.deepseek.com"))

    tasks = json.load(open(args.tasks))
    only = {t.strip() for t in args.only.split(",") if t.strip()}
    if only:
        tasks = [t for t in tasks if t["id"] in only]
        print(f"[only] running {len(tasks)} task(s): {[t['id'] for t in tasks]}", file=sys.stderr)
    os.makedirs(args.out_dir, exist_ok=True)
    modes = {"all": ["tools", "skill", "baseline"], "both": ["tools", "baseline"]}.get(
        args.mode, [args.mode])

    for mode in modes:
        rows = []
        for task in tasks:
            print(f"[{mode}] task {task['id']}: {task['prompt'][:70]}...", file=sys.stderr)
            r = run_task(client, args.model, task, mode)
            print(f"   -> success={r['success']} turns={r['turns']} "
                  f"tool_calls={r['tool_calls']} errors={r['tool_errors']} "
                  f"tok={r['prompt_tokens']+r['completion_tokens']} {r['wall_s']}s", file=sys.stderr)
            rows.append(r)

        suffix = ".only" if only else ""
        json.dump(rows, open(os.path.join(args.out_dir, f"{mode}{suffix}_transcripts.json"), "w"), indent=2)
        n = len(rows)
        agg = {
            "mode": mode, "model": args.model, "n_tasks": n,
            "success_rate": round(sum(r["success"] for r in rows) / n, 4) if n else 0,
            "avg_turns": round(sum(r["turns"] for r in rows) / n, 2) if n else 0,
            "avg_tool_calls": round(sum(r["tool_calls"] for r in rows) / n, 2) if n else 0,
            "total_tool_errors": sum(r["tool_errors"] for r in rows),
            "forced_answers": sum(1 for r in rows if r.get("forced_answer")),
            "avg_tokens": round(sum(r["prompt_tokens"] + r["completion_tokens"] for r in rows) / n, 0) if n else 0,
            "avg_wall_s": round(sum(r["wall_s"] for r in rows) / n, 2) if n else 0,
            "per_task": [{k: r[k] for k in ("task_id", "success", "turns", "tool_calls",
                                            "tool_errors", "forced_answer", "wall_s")} for r in rows],
        }
        json.dump(agg, open(os.path.join(args.out_dir, f"{mode}{suffix}.json"), "w"), indent=2)
        print(json.dumps({k: agg[k] for k in ("mode", "success_rate", "avg_turns",
              "avg_tool_calls", "total_tool_errors", "avg_tokens", "avg_wall_s")}, indent=2))

if __name__ == "__main__":
    main()
