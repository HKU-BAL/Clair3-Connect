#!/usr/bin/env bash

set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
source "${HERE}/env.sh"
source "${PAPER_EXP}/lib/common.sh"

N="${E4_REPS:-10}"
E4_PORT="${E4_PORT:-50075}"
APOE_TASKS="apoe_diplotype_hg00136,apoe_diplotype_hg00150"
OUT="${RESULTS}/e4_agent/repeats"; mkdir -p "$OUT"
export CLAIR3_CONNECT_SERVER="http://127.0.0.1:${E4_PORT}"
export CLAIR3_CONNECT_AGENT_WORK="${DATA_CACHE}/e4_agent_work"; mkdir -p "$CLAIR3_CONNECT_AGENT_WORK"

export REF_ONT REF_HIFI REF_ILMN CKPT_ONT CKPT_HIFI CKPT_ILMN CLAIR3_ROOT CLAIR3_CONNECT_BIN
export BCFTOOLS="${BCFTOOLS:-<CONDA>/envs/clair3-torch/bin/bcftools}"
export CLAIR_SKILL_DIR="${CLAIR_SKILL_DIR:-${DATA_CACHE}/Clair-skills}"

if [ -z "${DEEPSEEK_API_KEY:-}" ]; then
  echo "DEEPSEEK_API_KEY not set — export it at runtime and re-run." >&2
  exit 1
fi

if [ ! -f "${DATA_CACHE}/apoe/HG00136_apoe.bam" ] || [ ! -f "${DATA_CACHE}/apoe/HG00150_apoe.bam" ]; then
  bash "${HERE}/prep_apoe_data.sh" || log "APOE data prep failed"
fi
export APOE_TRUTH_JSON="${DATA_CACHE}/apoe/apoe_truth.json"
export APOE_BAM_HG00136="${DATA_CACHE}/apoe/HG00136_apoe.bam"

export SERVER_LOG="${LOGS}/e4_repeats_server.log"; : > "$SERVER_LOG"
start_server --bind "127.0.0.1:${E4_PORT}" --model-registry "$PAPER_REGISTRY" --device cpu --max-inflight 8
trap 'stop_server' EXIT
wait_server_ready "$SERVER_LOG" 180 || { log "server failed"; cat "$SERVER_LOG"; exit 1; }
log "server ready on $CLAIR3_CONNECT_SERVER (model ${DEEPSEEK_MODEL}, N=${N})"

source "${PAPER_EXP}/.venv-e4/bin/activate"

for i in $(seq 1 "$N"); do
  rep="$(printf 'rep%02d' "$i")"
  repdir="${OUT}/${rep}"; mkdir -p "$repdir"
  if [ -f "${repdir}/tools.only.json" ] && [ -f "${repdir}/skill.only.json" ] && [ -f "${repdir}/baseline.only.json" ]; then
    log "[$rep] already complete — skipping"
    continue
  fi
  log "[$rep] 3 arms × 2 APOE samples ..."
  if python "${HERE}/agent.py" --mode all --only "$APOE_TASKS" \
        --out-dir "$repdir" --model "$DEEPSEEK_MODEL" > "${LOGS}/e4_${rep}.log" 2>&1; then
    log "[$rep] done -> ${repdir}"
  else
    log "[$rep] FAILED (see ${LOGS}/e4_${rep}.log) — continuing"
  fi
done

log "all repeats finished -> ${OUT}"
