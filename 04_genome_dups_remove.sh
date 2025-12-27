#!/usr/bin/env bash
set -euo pipefail

RAW_DIR="./result/assembly/denovo/"
OUT_DIR="./result/clean/denovo"
SCRIPT="./tools/Remove_dups.py"
THREADS=32
COMMANDS="commandsCLEAN.cmd"

usage() {
cat <<EOF
Usage: $(basename "$0") [options]

Options:
  -i, --input DIR        Input FASTA directory (default: ./result/assembly/denovo)
  -o, --outdir DIR       Output directory (default: ./result/clean/denovo)
  -s, --script PATH      Remove_dups.py path (default: ./tools/Remove_dups.py)
  -t, --threads INT      Parallel jobs (default: 32)
  -c, --cmdfile FILE     Command list file (default: commandsCLEAN.cmd)
  -h, --help             Show this help

Example:
  $(basename "$0") \\
    -i ./result/assembly/denovo \\
    -o ./result/clean/denovo \\
    -s ./tools/Remove_dups.py \\
    -t 32

nohup:
  nohup $(basename "$0") -i ./result/assembly/denovo -o ./result/clean/denovo -t 32 \\
    > remove_dups.nohup.log 2>&1 &
EOF
}

#######################################
# parse args
#######################################
while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input)   RAW_DIR="$2"; shift 2 ;;
    -o|--outdir)  OUT_DIR="$2"; shift 2 ;;
    -s|--script)  SCRIPT="$2"; shift 2 ;;
    -t|--threads) THREADS="$2"; shift 2 ;;
    -c|--cmdfile) COMMANDS="$2"; shift 2 ;;
    -h|--help)    usage; exit 0 ;;
    *) echo "ERROR: Unknown option: $1"; usage; exit 1 ;;
  esac
done

#######################################
# checks
#######################################
if [[ ! -d "$RAW_DIR" ]]; then
  echo "❌ ERROR: input directory not found: $RAW_DIR"
  exit 1
fi

if [[ ! -f "$SCRIPT" ]]; then
  echo "❌ ERROR: script not found: $SCRIPT"
  exit 1
fi

if ! command -v parallel >/dev/null 2>&1; then
  echo "❌ ERROR: GNU parallel not found in PATH"
  exit 1
fi

mkdir -p "$OUT_DIR"
rm -f "$COMMANDS"

echo "=== Remove_dups batch started: $(date) ==="
echo "Input dir:   $RAW_DIR"
echo "Output dir:  $OUT_DIR"
echo "Script:      $SCRIPT"
echo "Threads:     $THREADS"
echo "Cmd file:    $COMMANDS"
echo "========================================="

#######################################
# generate command list
#######################################
shopt -s nullglob
FASTA_FILES=("$RAW_DIR"/*.fasta)
shopt -u nullglob

if [[ ${#FASTA_FILES[@]} -eq 0 ]]; then
  echo "❌ No .fasta files found in $RAW_DIR"
  exit 1
fi

for fasta in "${FASTA_FILES[@]}"; do
  name=$(basename "$fasta" .fasta)
  out="$OUT_DIR/$name"
  echo "python $SCRIPT -d $fasta -o $out" >> "$COMMANDS"
done

echo "Generated command list ($COMMANDS):"
cat "$COMMANDS"
echo "-----------------------------------"

#######################################
# run in parallel
#######################################
echo "Running Remove_dups using $THREADS threads..."
parallel -j "$THREADS" < "$COMMANDS"

echo "✅ All jobs finished at $(date)"
