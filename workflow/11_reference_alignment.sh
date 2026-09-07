#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/lib/common.sh"
validate_value_options "--output-dir --query --reference --threads -o -q -r -t" "$@"

# -------------------------
# Default values and variables
# -------------------------
REF=""
QRY=""
THREADS=32
MAXMATCH="--maxmatch"
NOSIMPLIFY="--nosimplify"
OUTPUT_DIR="."

# -------------------------
# help
# -------------------------
usage() {
  cat <<EOF
Usage:
  $(basename "$0") -r <ref.fa> -q <query.fasta> [options]

Required:
  -r, --reference   Reference fasta file (e.g., ref.fa)
  -q, --query       Query fasta file (e.g., query.fasta)

Options:
  -t, --threads     Number of threads to use (default: 32)
  -m, --maxmatch    Run with maxmatch (default: true)
  -n, --nosimplify  Run without simplify (default: true)
  -o, --output-dir  Directory to save .coords (default: current directory)
  -h, --help        Show this help message

Example:
  $(basename "$0") -r ref.fa -q query.fasta -t 16 -o /path/to/output

EOF
}

die() { echo "❌ $*" >&2; exit 1; }

# -------------------------
# Parse arguments
# -------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    -r|--reference)  REF="$2"; shift 2;;
    -q|--query)      QRY="$2"; shift 2;;
    -t|--threads)    THREADS="$2"; shift 2;;
    -m|--maxmatch)   MAXMATCH="--maxmatch"; shift;;
    -n|--nosimplify) NOSIMPLIFY="--nosimplify"; shift;;
    -o|--output-dir) OUTPUT_DIR="$2"; shift 2;;
    -h|--help)       usage; exit 0;;
    *) die "Unknown option: $1";;
  esac
done
[[ -z "${THREADS}" ]] || positive_int THREADS "${THREADS}"

# -------------------------
# Validate inputs
# -------------------------
[[ -f "$REF" ]] || die "Reference file not found: $REF"
[[ -f "$QRY" ]] || die "Query file not found: $QRY"
mkdir -p "$OUTPUT_DIR"

# -------------------------
# Processing
# -------------------------
# Get query basename and remove .fasta or .fa suffix
BASENAME=$(basename "$QRY")
PFX=${BASENAME%.fasta}
PFX=${PFX%.fa}
COORDS="$OUTPUT_DIR/$PFX.coords"
PFX="$OUTPUT_DIR/$PFX"

echo "[INFO] Reference : ${REF}"
echo "[INFO] Query     : ${QRY}"
echo "[INFO] Prefix    : ${PFX}"
echo "[INFO] Threads   : ${THREADS}"
echo "[INFO] Output dir: ${OUTPUT_DIR}"

# Run nucmer
nucmer \
  -t "$THREADS" \
  $MAXMATCH \
  $NOSIMPLIFY \
  -p "${PFX}" \
  "$REF" \
  "$QRY"

# Run delta-filter
delta-filter -m "${PFX}.delta" > "${PFX}.delta_filter"

# Generate coordinates
show-coords -c "${PFX}.delta_filter" > "${COORDS}"

# Remove temporary files
rm -f "${PFX}.delta" "${PFX}.delta_filter"

echo "[DONE] Outputs:"
echo "  ${COORDS}"
