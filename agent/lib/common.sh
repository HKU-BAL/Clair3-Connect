#!/usr/bin/env bash

set -o pipefail

log()  { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*" >&2; }
die()  { printf '[FATAL] %s\n' "$*" >&2; exit 1; }

wait_server_ready() {
  local logf="$1" timeout_s="${2:-120}"
  timeout "$timeout_s" grep -q -m1 -E 'listening on' \
      < <(tail -n +1 -f "$logf" 2>/dev/null)
}

wait_port() {
  local host="$1" port="$2" timeout_s="${3:-120}"
  timeout "$timeout_s" bash -c \
    'until (exec 3<>/dev/tcp/'"$0"'/'"$1"') 2>/dev/null; do :; done' "$host" "$port"
}

start_server() {
  local logf="${SERVER_LOG:-${LOGS}/server_$$.log}"
  env -u LD_LIBRARY_PATH "$CLAIR3_CONNECT_SERVER_BIN" "$@" --non-interactive \
      >"$logf" 2>&1 &
  SERVER_PID=$!
  log "server pid=$SERVER_PID log=$logf args: $*"
}

stop_server() {
  local pid="${1:-$SERVER_PID}"
  [ -n "$pid" ] && kill "$pid" 2>/dev/null
  wait "$pid" 2>/dev/null
}

now_ns() { date +%s%N; }

stamp() {
  local exp="$1" out="$2"; shift 2
  activate_torch 2>/dev/null
  python "${PAPER_EXP}/lib/provenance.py" --experiment "$exp" --out "$out" \
      --extra "$@" 2>/dev/null || true
}
