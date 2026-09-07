#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/lib/common.sh"
validate_value_options "--assembly-dir --kmer --memory --outdir --pattern --r1-regex --r2-regex --reads-dir --suffix --threads -a -k -m -o -r -t" "$@"

# ==========================================================
# Batch Merqury QV estimation (Illumina PE)
# ==========================================================

# -------------------------
# defaults
# -------------------------
READS_DIR=""
ASSEMBLY_DIR=""
WORKDIR=""

KMER=17
THREADS=4
MEMORY="8g"

ASM_GLOB="*.nextpolish.fasta"
ASM_SUFFIX=".nextpolish.fasta"

R1_REGEX='[_.-]1\.fq\.gz$'
R2_REGEX='[_.-]2\.fq\.gz$'

DRY_RUN=false
SKIP_EXISTING=true

# -------------------------
# help
# -------------------------
usage() {
cat <<EOF
Usage:
  $(basename "$0") -a ASSEMBLY_DIR -r READS_DIR -o WORKDIR [options]

Required:
  -a, --assembly-dir DIR     Assembly fasta directory
  -r, --reads-dir DIR        Illumina reads directory
  -o, --outdir DIR           Output working directory

Options:
  -k, --kmer INT             k-mer size (default: 17)
  -t, --threads INT          Threads for meryl (default: 4)
  -m, --memory STR           Memory for meryl (default: 8g)

  --pattern GLOB             Assembly glob (default: *.nextpolish.fasta)
  --suffix STR               Assembly suffix to strip for PREFIX

  --r1-regex REGEX           Regex for R1 (default: [_.-]1\\.fq\\.gz$)
  --r2-regex REGEX           Regex for R2 (default: [_.-]2\\.fq\\.gz$)

  --dry-run                  Print commands only, do not execute
  --no-skip-existing         Re-run even if outputs exist
  -h, --help                 Show this help

Naming logic:
  PREFIX = assembly filename stripped by --suffix
  Reads detected by:
    .*/PREFIX.*<R1_REGEX>
    .*/PREFIX.*<R2_REGEX>

Example:
  $(basename "$0") \\
    -a ./assembly_clean \\
    -r ./illumina \\
    -o ./merqury_qv \\
    -k 17 -t 8 -m 16g

EOF
}

die() { echo "❌ $*" >&2; exit 1; }
run() {
  if [[ "$DRY_RUN" == true ]]; then
    printf "[DRY-RUN] "; printf "%q " "$@"; printf "\n"
  else
    "$@"
  fi
}

# -------------------------
# arg parsing
# -------------------------
if [[ $# -eq 0 ]]; then usage; exit 1; fi

while [[ $# -gt 0 ]]; do
  case "$1" in
    -a|--assembly-dir) ASSEMBLY_DIR="$2"; shift 2;;
    -r|--reads-dir)    READS_DIR="$2"; shift 2;;
    -o|--outdir)       WORKDIR="$2"; shift 2;;

    -k|--kmer)         KMER="$2"; shift 2;;
    -t|--threads)      THREADS="$2"; shift 2;;
    -m|--memory)       MEMORY="$2"; shift 2;;

    --pattern)         ASM_GLOB="$2"; shift 2;;
    --suffix)          ASM_SUFFIX="$2"; shift 2;;
    --r1-regex)        R1_REGEX="$2"; shift 2;;
    --r2-regex)        R2_REGEX="$2"; shift 2;;

    --dry-run)         DRY_RUN=true; shift;;
    --no-skip-existing) SKIP_EXISTING=false; shift;;

    -h|--help) usage; exit 0;;
    *) die "Unknown option: $1";;
  esac
done
[[ -z "${THREADS}" ]] || positive_int THREADS "${THREADS}"
[[ -z "${KMER}" ]] || positive_int KMER "${KMER}"

# -------------------------
# validate
# -------------------------
[[ -d "$ASSEMBLY_DIR" ]] || die "Assembly dir not found"
[[ -d "$READS_DIR" ]]    || die "Reads dir not found"
[[ -n "$WORKDIR" ]]      || die "Missing --outdir"

mkdir -p "$WORKDIR"

ASSEMBLY_DIR=$(realpath "$ASSEMBLY_DIR")
READS_DIR=$(realpath "$READS_DIR")
WORKDIR=$(realpath "$WORKDIR")

echo "=============================="
echo "Merqury QV configuration"
echo "  assembly_dir : $ASSEMBLY_DIR"
echo "  reads_dir    : $READS_DIR"
echo "  outdir       : $WORKDIR"
echo "  k-mer        : $KMER"
echo "  threads      : $THREADS"
echo "  memory       : $MEMORY"
echo "  asm_glob     : $ASM_GLOB"
echo "  skip_existing: $SKIP_EXISTING"
echo "=============================="

# -------------------------
# main loop
# -------------------------
shopt -s nullglob
assemblies=( "$ASSEMBLY_DIR"/$ASM_GLOB )
(( ${#assemblies[@]} > 0 )) || die "No assemblies matched"

for fasta in "${assemblies[@]}"; do
  base=$(basename "$fasta")
  PREFIX="${base%$ASM_SUFFIX}"

  echo "=== Processing $PREFIX ==="

  outdir="$WORKDIR/${PREFIX}.merquryQV"
  if [[ "$SKIP_EXISTING" == true && -d "$outdir" && -f "$outdir/${PREFIX}.out.qv" ]]; then
    echo "[SKIP] $outdir exists and ${PREFIX}.out.qv found"
    continue
  fi
  mkdir -p "$outdir"
  cd "$outdir"

  R1=$(find_unique_read "$READS_DIR" "$PREFIX" "$R1_REGEX")
  R2=$(find_unique_read "$READS_DIR" "$PREFIX" "$R2_REGEX")

  if [[ -z "$R1" || -z "$R2" ]]; then
    echo "⚠️  WARNING: R1/R2 not found for $PREFIX → skipping"
    continue
  fi

  echo "[INFO] R1: $R1"
  echo "[INFO] R2: $R2"

  if [[ ! -d "${PREFIX}.meryl" || "$SKIP_EXISTING" == false ]]; then
    run meryl k=$KMER count memory=$MEMORY threads=$THREADS \
        output ${PREFIX}.meryl "$R1" "$R2"
  else
    echo "[SKIP] ${PREFIX}.meryl exists"
  fi

  if [[ ! -f "${PREFIX}.out.qv" || "$SKIP_EXISTING" == false ]]; then
    run merqury.sh "${PREFIX}.meryl" "$fasta" "${PREFIX}.out"
  else
    echo "[SKIP] ${PREFIX}.out.qv exists"
  fi

  echo "=== Finished $PREFIX ==="
  echo
done

# -------------------------
# summarize
# -------------------------
[[ "$DRY_RUN" == false ]] || exit 0
cd "$WORKDIR"
echo -e "Sample\tQV" > Collapsed.qv

for f in */*.out.qv; do
  sample=$(basename "$(dirname "$f")" .merquryQV)
  qv=$(awk '{print $4}' "$f" | head -n1)
  echo -e "${sample}\t${qv}" >> Collapsed.qv
done

echo "✅ All done! QV summary: $WORKDIR/Collapsed.qv"
