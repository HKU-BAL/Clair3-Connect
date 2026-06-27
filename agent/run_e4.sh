#!/usr/bin/env bash

set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
source "${HERE}/env.sh"
source "${PAPER_EXP}/lib/common.sh"

E4_PORT="${E4_PORT:-50074}"
OUT="${RESULTS}/e4_agent"; mkdir -p "$OUT"
export CLAIR3_CONNECT_SERVER="http://127.0.0.1:${E4_PORT}"
export CLAIR3_CONNECT_AGENT_WORK="${DATA_CACHE}/e4_agent_work"; mkdir -p "$CLAIR3_CONNECT_AGENT_WORK"

export REF_ONT REF_HIFI REF_ILMN CKPT_ONT CKPT_HIFI CKPT_ILMN CLAIR3_ROOT CLAIR3_CONNECT_BIN
export BCFTOOLS="${BCFTOOLS:-<CONDA>/envs/clair3-torch/bin/bcftools}"

export CLAIR_SKILL_DIR="${CLAIR_SKILL_DIR:-${DATA_CACHE}/Clair-skills}"
if [ ! -f "${CLAIR_SKILL_DIR}/SKILL.md" ]; then
  log "cloning HKU-BAL/Clair-skills -> ${CLAIR_SKILL_DIR}"
  git clone --depth 1 https://github.com/HKU-BAL/Clair-skills "${CLAIR_SKILL_DIR}" >/dev/null 2>&1 || log "clone failed (set CLAIR_SKILL_DIR)"
fi

bash "${HERE}/prep_apoe_data.sh" || log "APOE data prep failed (APOE task will error)"
export APOE_TRUTH_JSON="${DATA_CACHE}/apoe/apoe_truth.json"

export SERVER_LOG="${LOGS}/e4_server.log"; : > "$SERVER_LOG"
start_server --bind "127.0.0.1:${E4_PORT}" --model-registry "$PAPER_REGISTRY" --device cpu --max-inflight 8
trap 'stop_server' EXIT
wait_server_ready "$SERVER_LOG" 180 || { log "server failed"; cat "$SERVER_LOG"; exit 1; }
log "server ready on $CLAIR3_CONNECT_SERVER"

source "${PAPER_EXP}/.venv-e4/bin/activate"

if [ -z "${DEEPSEEK_API_KEY:-}" ]; then
  log "DEEPSEEK_API_KEY not set — set it (e.g.  ! export DEEPSEEK_API_KEY=sk-...) and re-run."
  exit 0
fi
log "running DeepSeek agent (tools + skill + baseline) with model ${DEEPSEEK_MODEL}"
python "${HERE}/agent.py" --mode "${E4_MODE:-all}" --out-dir "$OUT" --model "$DEEPSEEK_MODEL"
stamp e4_agent "${OUT}/manifest.json" "model=${DEEPSEEK_MODEL}" "server=cpu"
log "E4 done: $(ls "$OUT"/*.json 2>/dev/null | tr '\n' ' ')"
