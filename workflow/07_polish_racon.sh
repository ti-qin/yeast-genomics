#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/lib/common.sh"
validate_value_options "--assembly-dir --minimap2 --outdir --pattern --racon --reads-dir --reads-suffix --rounds --suffix --threads -a -n -o -r -t" "$@"

# ==========================================================
# Batch multi-round Racon polishing (ONT)
# ==========================================================

# -------------------------
# defaults
# -------------------------
ASSEMBLY_DIR=""
READS_DIR=""
OUT_DIR=""

THREADS=96
ROUNDS=3

ASM_GLOB="*.nextdenovo.NR.nuclear.fasta"
ASM_SUFFIX=".nextdenovo.NR.nuclear.fasta"
READS_SUFFIX=".fastq.gz"

MINIMAP2="minimap2"
RACON="racon"

KEEP_WORKDIR=false

# -------------------------
# help
# -------------------------
usage() {
cat <<EOF
Usage:
  $(basename "$0") -a ASSEMBLY_DIR -r READS_DIR -o OUT_DIR [options]

Required:
  -a, --assembly-dir DIR     Directory containing assemblies
  -r, --reads-dir DIR        Directory containing ONT reads
  -o, --outdir DIR           Output directory

Options:
  -t, --threads INT          Threads (default: 96)
  -n, --rounds INT           Number of Racon rounds (default: 3)

  --pattern GLOB             Assembly glob (default: *.nextdenovo.NR.nuclear.fasta)
  --suffix STR               Assembly suffix to strip for PREFIX
                             (default: .nextdenovo.NR.nuclear.fasta)

  --reads-suffix STR         Reads filename suffix (default: .fastq.gz)

  --minimap2 PATH            minimap2 executable (default: minimap2)
  --racon PATH               racon executable (default: racon)

  --keep-workdir             Keep per-sample mapping directory
  -h, --help                 Show this help

Naming rule:
  PREFIX is inferred from assembly filename by stripping --suffix.
  Reads must be:
    <READS_DIR>/<PREFIX><reads-suffix>

Example:
  $(basename "$0") \\
    -a ./assemblies \\
    -r ./reads \\
    -o ./racon_polished \\
    -t 96 -n 3

EOF
}

die() { echo "❌ $*" >&2; exit 1; }

# -------------------------
# arg parsing
# -------------------------
if [[ $# -eq 0 ]]; then usage; exit 1; fi

while [[ $# -gt 0 ]]; do
  case "$1" in
    -a|--assembly-dir) ASSEMBLY_DIR="$2"; shift 2;;
    -r|--reads-dir)    READS_DIR="$2"; shift 2;;
    -o|--outdir)       OUT_DIR="$2"; shift 2;;

    -t|--threads)      THREADS="$2"; shift 2;;
    -n|--rounds)       ROUNDS="$2"; shift 2;;

    --pattern)         ASM_GLOB="$2"; shift 2;;
    --suffix)          ASM_SUFFIX="$2"; shift 2;;
    --reads-suffix)    READS_SUFFIX="$2"; shift 2;;

    --minimap2)        MINIMAP2="$2"; shift 2;;
    --racon)           RACON="$2"; shift 2;;

    --keep-workdir)    KEEP_WORKDIR=true; shift;;

    -h|--help) usage; exit 0;;
    *) die "Unknown option: $1";;
  esac
done
[[ -z "${THREADS}" ]] || positive_int THREADS "${THREADS}"
[[ -z "${ROUNDS}" ]] || positive_int ROUNDS "${ROUNDS}"

# -------------------------
# validate
# -------------------------
[[ -n "$ASSEMBLY_DIR" ]] || die "Missing --assembly-dir"
[[ -n "$READS_DIR" ]]    || die "Missing --reads-dir"
[[ -n "$OUT_DIR" ]]      || die "Missing --outdir"

[[ -d "$ASSEMBLY_DIR" ]] || die "Assembly dir not found"
[[ -d "$READS_DIR" ]]    || die "Reads dir not found"

command -v "$MINIMAP2" >/dev/null || die "minimap2 not found: $MINIMAP2"
command -v "$RACON" >/dev/null || die "racon not found: $RACON"

mkdir -p "$OUT_DIR"

ASSEMBLY_DIR="$(realpath "$ASSEMBLY_DIR")"
READS_DIR="$(realpath "$READS_DIR")"
OUT_DIR="$(realpath "$OUT_DIR")"

echo "=============================="
echo "Racon polishing configuration"
echo "  assembly_dir : $ASSEMBLY_DIR"
echo "  reads_dir    : $READS_DIR"
echo "  outdir       : $OUT_DIR"
echo "  threads      : $THREADS"
echo "  rounds       : $ROUNDS"
echo "  asm_glob     : $ASM_GLOB"
echo "  asm_suffix   : $ASM_SUFFIX"
echo "  reads_suffix : $READS_SUFFIX"
echo "  keep_workdir : $KEEP_WORKDIR"
echo "=============================="

# -------------------------
# main loop
# -------------------------
shopt -s nullglob
assemblies=( "$ASSEMBLY_DIR"/$ASM_GLOB )
(( ${#assemblies[@]} > 0 )) || die "No assemblies matched"

for ASSEMBLY in "${assemblies[@]}"; do
  base="$(basename "$ASSEMBLY")"
  PREFIX="${base%$ASM_SUFFIX}"

  READS="$READS_DIR/${PREFIX}${READS_SUFFIX}"
  [[ -f "$READS" ]] || {
    die "Reads not found for $PREFIX: $READS"
  }

  echo "========== [$PREFIX] Racon polishing =========="

  WORKDIR="$OUT_DIR/${PREFIX}_NanoporeMapping"
  mkdir -p "$WORKDIR"

  INPUT="$ASSEMBLY"

  for ((ROUND=1; ROUND<=ROUNDS; ROUND++)); do
    echo "[$(date)] Round $ROUND mapping"
    PAF="$WORKDIR/${PREFIX}_r${ROUND}.paf"

    "$MINIMAP2" -x map-ont -t "$THREADS" "$INPUT" "$READS" > "$PAF"

    echo "[$(date)] Round $ROUND polishing"
    OUTPUT="$WORKDIR/${PREFIX}_r${ROUND}.fasta"

    "$RACON" -t "$THREADS" "$READS" "$PAF" "$INPUT" > "$OUTPUT"

    INPUT="$OUTPUT"
  done

  FINAL_OUT="$OUT_DIR/${PREFIX}.recon${ROUNDS}.fasta"
  cp "$INPUT" "$FINAL_OUT"

  echo "[$(date)] [$PREFIX] Finished → $FINAL_OUT"

  if [[ "$KEEP_WORKDIR" == false ]]; then
    rm -rf "$WORKDIR"
  fi
done

echo "🎉 All Racon polishing jobs completed!"
