#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/lib/common.sh"
validate_value_options "--indir --out -i -o" "$@"

# ---------------------------
# defaults
# ---------------------------
INDIR=""
OUTFILE="assembly_summary.tsv"

# ---------------------------
# usage
# ---------------------------
usage() {
cat <<EOF
Usage: $(basename "$0") -i DIR [options]

Required:
  -i, --indir DIR        directory containing *.fasta files

Optional:
  -o, --out FILE         output summary file (default: assembly_summary.tsv)
  -h, --help             show this help message

Example:
  $(basename "$0") -i ./result/clean/denovo -o assembly_summary_denovo.tsv

EOF
}

# ---------------------------
# parse args
# ---------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--indir)
      INDIR="$2"
      shift 2
      ;;
    -o|--out)
      OUTFILE="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "❌ ERROR: Unknown option: $1"
      usage
      exit 1
      ;;
  esac
done

# ---------------------------
# checks
# ---------------------------
if [[ -z "$INDIR" ]]; then
  echo "❌ ERROR: --indir is required"
  usage
  exit 1
fi

if [[ ! -d "$INDIR" ]]; then
  echo "❌ ERROR: Directory not found: $INDIR"
  exit 1
fi

if ! command -v assembly-stats >/dev/null 2>&1; then
  echo "❌ ERROR: assembly-stats not found in PATH"
  exit 1
fi

# ---------------------------
# header
shopt -s nullglob
FASTAS=( "$INDIR"/*.fasta )
(( ${#FASTAS[@]} > 0 )) || { echo "ERROR: No *.fasta inputs" >&2; exit 1; }
mkdir -p "$(dirname "$OUTFILE")"
# ---------------------------
echo -e "Sample\tSum_bp\tNum_seqs\tAverage\tLargest\tN50\tN90" > "$OUTFILE"

# ---------------------------
# main loop
# ---------------------------
shopt -s nullglob
for f in "$INDIR"/*.fasta; do
  [[ -f "$f" ]] || continue

  name=$(basename "$f" .fasta)
  stats=$(assembly-stats "$f")

  sum=$(echo "$stats" | grep "^sum" | awk -F'[ =,]+' '{print $2}')
  nseq=$(echo "$stats" | grep "^sum" | awk -F'[ =,]+' '{print $4}')
  ave=$(echo "$stats" | grep "^sum" | awk -F'[ =,]+' '{print $6}')
  largest=$(echo "$stats" | grep "^sum" | awk -F'[ =,]+' '{print $8}')
  n50=$(echo "$stats" | grep "^N50" | awk -F'[ =,]+' '{print $2}')
  n90=$(echo "$stats" | grep "^N90" | awk -F'[ =,]+' '{print $2}')

  echo -e "${name}\t${sum}\t${nseq}\t${ave}\t${largest}\t${n50}\t${n90}" >> "$OUTFILE"
done
shopt -u nullglob

echo "✅ Summary saved to $OUTFILE"
