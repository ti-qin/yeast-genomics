#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  15_mash_distance.sh <query_fasta_dir> <ref_db_msh> <out_dir>

Example:
  bash 15_mash_distance.sh ../result/final_asm ref332db.msh ./mash_out

Outputs:
  out_dir/
    sketches/            (generated *.msh)
    G9148_vs_ref332db.res
    G1234_vs_ref332db.res
    ...
EOF
}

if [[ "${1:-}" == -h || "${1:-}" == --help ]]; then usage; exit 0; fi

if [[ $# -ne 3 ]]; then
  usage
  exit 1
fi

QUERY_DIR="$1"
REF_DB_MSH="$2"
OUT_DIR="$3"

# Params (you can change if needed)
K=21
S=10000
THREADS=${THREADS:-1}   # Mash sketch is mostly single-thread; keep for future use

# Basic checks
if [[ ! -d "$QUERY_DIR" ]]; then
  echo "[ERROR] query_fasta_dir not found: $QUERY_DIR" >&2
  exit 2
fi
if [[ ! -f "$REF_DB_MSH" ]]; then
  echo "[ERROR] ref_db_msh not found: $REF_DB_MSH" >&2
  exit 2
fi
command -v mash >/dev/null 2>&1 || { echo "[ERROR] mash not in PATH" >&2; exit 2; }

mkdir -p "$OUT_DIR"
SKETCH_DIR="${OUT_DIR%/}/sketches"
mkdir -p "$SKETCH_DIR"

# Make a safe tag from ref db filename (for output naming)
REF_TAG="$(basename "$REF_DB_MSH")"
REF_TAG="${REF_TAG%.msh}"

shopt -s nullglob

# Collect fasta-like files
FILES=( "$QUERY_DIR"/*.fa "$QUERY_DIR"/*.fasta "$QUERY_DIR"/*.fna \
        "$QUERY_DIR"/*.fa.gz "$QUERY_DIR"/*.fasta.gz "$QUERY_DIR"/*.fna.gz )

if [[ ${#FILES[@]} -eq 0 ]]; then
  echo "[ERROR] No fasta files found in: $QUERY_DIR" >&2
  exit 3
fi

echo "[INFO] Query dir : $QUERY_DIR"
echo "[INFO] Ref db    : $REF_DB_MSH"
echo "[INFO] Out dir   : $OUT_DIR"
echo "[INFO] Sketch dir: $SKETCH_DIR"
echo "[INFO] Found ${#FILES[@]} fasta files"

declare -A SEEN_SAMPLES=()
for f in "${FILES[@]}"; do
  base="$(basename "$f")"

  # Derive sample name:
  # - If filename like G9148.nextpolish.fasta -> sample = G9148
  # - Else use filename without extension(s)
  sample="${base%%.*}"
  if [[ -z "$sample" || "$sample" == "$base" ]]; then
    sample="$base"
    sample="${sample%.gz}"
    sample="${sample%.fasta}"
    sample="${sample%.fa}"
    sample="${sample%.fna}"
  fi

  [[ -z "${SEEN_SAMPLES[$sample]:-}" ]] || { echo "ERROR: Multiple assemblies map to sample $sample" >&2; exit 1; }
  SEEN_SAMPLES[$sample]=1
  sketch_prefix="${SKETCH_DIR}/${sample}"
  sketch_file="${sketch_prefix}.msh"
  out_file="${OUT_DIR%/}/${sample}_vs_${REF_TAG}.res"

  echo "==> [$sample] input: $f"
  # If gz, mash can read via process substitution
  if [[ "$f" == *.gz ]]; then
    mash sketch -o "$sketch_prefix" -k "$K" -s "$S" <(gunzip -c "$f")
  else
    mash sketch -o "$sketch_prefix" -k "$K" -s "$S" "$f"
  fi

  mash dist "$sketch_file" "$REF_DB_MSH" > "$out_file"
  echo "    wrote: $out_file"
done

echo "[DONE] All finished. Results in: $OUT_DIR"
