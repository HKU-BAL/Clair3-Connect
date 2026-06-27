#!/usr/bin/env bash

set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
source "${HERE}/env.sh"

APOE_DIR="${APOE_DIR:-${DATA_CACHE}/apoe}"; mkdir -p "$APOE_DIR"

APOE_REGION="${APOE_REGION:-chr19:44903000-44913000}"
SAMTOOLS="${SAMTOOLS:-<CONDA>/envs/clair3-torch/bin/samtools}"

APOE_REF="${APOE_REF:-$REF_ILMN}"

export REF_CACHE="${APOE_DIR}/ref_cache/%2s/%2s/%s"

declare -A APOE_SAMPLES=(
  [HG00136]="https://ftp.sra.ebi.ac.uk/vol1/run/ERR324/ERR3240142/HG00136.final.cram;ilmn;e4/e4;1/1;0/0"
  [HG00150]="https://ftp.sra.ebi.ac.uk/vol1/run/ERR324/ERR3240154/HG00150.final.cram;ilmn;e2/e4;0/1;0/1"
)

TRUTH_JSON="${APOE_DIR}/apoe_truth.json"
echo "{" > "$TRUTH_JSON"
first=1
for s in "${!APOE_SAMPLES[@]}"; do
  IFS=';' read -r cram plat diplo gt1 gt2 <<< "${APOE_SAMPLES[$s]}"
  bam="${APOE_DIR}/${s}_apoe.bam"
  if [ ! -s "$bam" ] || [ ! -s "${bam}.bai" ]; then
    echo "[prep_apoe] streaming $s APOE region from 1000G CRAM ..."
    if timeout 300 "$SAMTOOLS" view -T "$APOE_REF" -b "$cram" "$APOE_REGION" -o "$bam" 2>"${APOE_DIR}/${s}.stream.err"; then
      "$SAMTOOLS" index "$bam"
    else
      echo "[prep_apoe] WARN: stream failed for $s (see ${APOE_DIR}/${s}.stream.err)"; continue
    fi
  fi
  n=$("$SAMTOOLS" view -c "$bam" 2>/dev/null || echo 0)
  echo "[prep_apoe] $s: $n reads, truth=$diplo  -> $bam"
  [ $first -eq 0 ] && echo "," >> "$TRUTH_JSON"; first=0
  printf '  "%s": {"platform": "%s", "diplotype": "%s", "rs429358_gt": "%s", "rs7412_gt": "%s", "bam": "%s"}' \
    "$s" "$plat" "$diplo" "$gt1" "$gt2" "$bam" >> "$TRUTH_JSON"
done
echo "" >> "$TRUTH_JSON"; echo "}" >> "$TRUTH_JSON"
echo "[prep_apoe] wrote $TRUTH_JSON"
cat "$TRUTH_JSON"
