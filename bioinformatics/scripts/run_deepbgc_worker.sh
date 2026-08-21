#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 GENOME_FASTA [OUTPUT_DIR]" >&2
  exit 2
fi

genome="$1"
output_root="${2:-deepbgc_results}"
name="$(basename "$genome")"
name="${name%.*}"

mkdir -p "$output_root"
echo "[$(date +"%T")] Processing $name"

if deepbgc pipeline "$genome" --output-dir "$output_root/$name" > "$output_root/${name}.log" 2>&1; then
  echo "[$(date +"%T")] Done $name"
else
  echo "[$(date +"%T")] FAILED: $name (see $output_root/${name}.log)" >&2
  exit 1
fi
