#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/lib/common.sh"
validate_value_options "--busco-lib --input-dir --jobs --lineage --outdir --suffix --threads -L -i -j -l -o -t" "$@"

# ==========================================================
# Batch Compleasm (genome mode)
# ==========================================================

# -------------------------
# defaults
# -------------------------
ASSEMBLY_DIR=""
WORKDIR=""

BUSCO_LINEAGE="saccharomycetes"
BUSCO_LIB=""
BUSCO_THREADS=8
PARALLEL_JOBS=4

PROT_SUFFIX=".fasta"
KEEP_WORKDIR=true

# -------------------------
# help
# -------------------------
usage() {
cat <<EOF
Usage:
  $(basename "$0") -i ASSEMBLY_DIR -o WORKDIR -L BUSCO_LIB [options]

Required:
  -i, --input-dir DIR        Directory containing genome FASTA (*.fasta)
  -o, --outdir DIR           Working/output directory
  -L, --busco-lib DIR        BUSCO lineage library directory

Options:
  -l, --lineage STR          BUSCO lineage name (default: saccharomycetes)
  -t, --threads INT          Threads per BUSCO run (default: 8)
  -j, --jobs INT             Parallel samples to run (default: 4)

  --suffix STR               Genome FASTA suffix (default: .fasta)
  --no-keep-workdir          Remove per-sample output after summary
  -h, --help                 Show this help

Behavior:
  - Runs: compleasm run
  - Skips samples if <sample>.miniBusco/summary.txt exists
  - Summarizes S/D/F/M and S+D score

Example:
  $(basename "$0") \\
    -i ./assemblies \\
    -o ./busco_res \\
    -L /path/to/BUSCO_lib \\
    -l saccharomycetes \\
    -t 8 -j 4

EOF
}

die() { echo "❌ $*" >&2; exit 1; }

# -------------------------
# arg parsing
# -------------------------
if [[ $# -eq 0 ]]; then usage; exit 1; fi

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input-dir)  ASSEMBLY_DIR="$2"; shift 2;;
    -o|--outdir)     WORKDIR="$2"; shift 2;;
    -L|--busco-lib)  BUSCO_LIB="$2"; shift 2;;

    -l|--lineage)    BUSCO_LINEAGE="$2"; shift 2;;
    -t|--threads)    BUSCO_THREADS="$2"; shift 2;;
    -j|--jobs)       PARALLEL_JOBS="$2"; shift 2;;

    --suffix)        PROT_SUFFIX="$2"; shift 2;;
    --no-keep-workdir) KEEP_WORKDIR=false; shift;;

    -h|--help) usage; exit 0;;
    *) die "Unknown option: $1";;
  esac
done
[[ -z "${PARALLEL_JOBS}" ]] || positive_int PARALLEL_JOBS "${PARALLEL_JOBS}"
[[ -z "${BUSCO_THREADS}" ]] || positive_int BUSCO_THREADS "${BUSCO_THREADS}"

# -------------------------
# validate
# -------------------------
[[ -d "$ASSEMBLY_DIR" ]] || die "Input dir not found"
[[ -n "$WORKDIR" ]]      || die "Missing --outdir"
[[ -d "$BUSCO_LIB" ]]    || die "BUSCO lib not found"

command -v compleasm >/dev/null || die "compleasm not found in PATH"

mkdir -p "$WORKDIR"

ASSEMBLY_DIR=$(realpath "$ASSEMBLY_DIR")
WORKDIR=$(realpath "$WORKDIR")
BUSCO_LIB=$(realpath "$BUSCO_LIB")

echo "=============================="
echo "MiniBUSCO / Compleasm config"
echo "  input_dir   : $ASSEMBLY_DIR"
echo "  outdir      : $WORKDIR"
echo "  lineage     : $BUSCO_LINEAGE"
echo "  busco_lib   : $BUSCO_LIB"
echo "  threads     : $BUSCO_THREADS"
echo "  parallel    : $PARALLEL_JOBS"
echo "  suffix      : $PROT_SUFFIX"
echo "=============================="

# -------------------------
# worker function
# -------------------------
run_minibusco() {
  local ASSEMBLY="$1"
  local SAMPLE
  SAMPLE=$(basename "$ASSEMBLY" "$PROT_SUFFIX")

  local OUTDIR="$WORKDIR/${SAMPLE}.miniBusco"
  local SUMMARY="$OUTDIR/summary.txt"

  if [[ -f "$SUMMARY" ]]; then
    echo "[SKIP] $SAMPLE already done"
    return
  fi

  echo "[RUN ] $SAMPLE"
  compleasm run \
    -a "$ASSEMBLY" \
    -o "$OUTDIR" \
    -L "$BUSCO_LIB" \
    -l "$BUSCO_LINEAGE" \
    -t "$BUSCO_THREADS"
}

export -f run_minibusco
export WORKDIR BUSCO_LIB BUSCO_LINEAGE BUSCO_THREADS PROT_SUFFIX

# -------------------------
# run in parallel
shopt -s nullglob
INPUTS=( "$ASSEMBLY_DIR"/*"$PROT_SUFFIX" )
(( ${#INPUTS[@]} > 0 )) || die "No matching FASTA files"
# -------------------------
find "$ASSEMBLY_DIR" -maxdepth 1 -type f -name "*${PROT_SUFFIX}" -print0 | \
  xargs -0 -r -P "$PARALLEL_JOBS" -I {} bash -c 'run_minibusco "$@"' _ {}

# -------------------------
# summarize
# -------------------------
cd "$WORKDIR"
echo -e "Sample\tBUSCO_Score\tS\tD\tF\tM" > BuscoScores.tsv

for F in *.miniBusco/summary.txt; do
  SAMPLE=$(basename "$(dirname "$F")" .miniBusco)

  S=$(grep "^S" "$F" | cut -d ":" -f2 | cut -d "%" -f1 | xargs)
  D=$(grep "^D" "$F" | cut                                                -d ":" -f2 | cut -d "%" -f1 | xargs)
  Fg=$(grep "^F" "$F" | cut -d ":" -f2 | cut -d "%" -f1 | xargs)
  M=$(grep "^M" "$F" | cut -d ":" -f2 | cut -d "%" -f1 | xargs)

  SCORE=$(awk "BEGIN{print $S + $D}")
  echo -e "${SAMPLE}\t${SCORE}\t${S}\t${D}\t${Fg}\t${M}"
done >> BuscoScores.tsv

echo "✅ BUSCO summary written to: $WORKDIR/BuscoScores.tsv"

if [[ "$KEEP_WORKDIR" == false ]]; then
  for summary in "$WORKDIR"/*.miniBusco/summary.txt; do
    rm -rf -- "$(dirname "$summary")"
  done
fi
