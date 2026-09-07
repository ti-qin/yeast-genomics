#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/lib/common.sh"
validate_value_options "--genome-dir --image --outdir --rnaseq-dir --suffix --threads -g -i -o -r -t" "$@"

# ==========================================================
# Batch BRaKER3 (RNA-seq) with Singularity
# ==========================================================

# -------------------------
# defaults
# -------------------------
FASTA_DIR=""
RNA_SEQ_DIR=""
OUTDIR="./braker_res"

SINGULARITY_IMAGE="braker3.sif"
THREADS=128

GENOME_SUFFIX=".masked.fasta"

BIND_RNA=true
NO_CLEANUP=true
FUNGUS=true
AB_INITIO=true
GFF3=true

# -------------------------
# help
# -------------------------
usage() {
cat <<EOF
Usage:
  $(basename "$0") -g GENOME_DIR -r RNA_SEQ_DIR -o OUTDIR [options]

Required:
  -g, --genome-dir DIR       Directory with genome FASTA files
  -r, --rnaseq-dir DIR       RNA-seq directory (for BRaKER3)
  -o, --outdir DIR           Output base directory

Options:
  -i, --image FILE           Singularity image (default: braker3.sif)
  -t, --threads INT          Threads (default: 128)

  --suffix STR               Genome FASTA suffix (default: .masked.fasta)

  --no-rna-bind              Do not bind RNA-seq dir explicitly
  --cleanup                  Allow BRaKER cleanup (default: --nocleanup)

  --no-fungus                Disable --fungus
  --no-ab-initio             Disable --AUGUSTUS_ab_initio
  --no-gff3                  Disable --gff3

  -h, --help                 Show this help

Example:
  $(basename "$0") \\
    -g ./braker_res/genome_fa \\
    -r ./rawdata/rnaseq \\
    -o ./braker_res \\
    -i braker3.sif \\
    -t 32

EOF
}

die() { echo "❌ $*" >&2; exit 1; }

# -------------------------
# arg parsing
# -------------------------
if [[ $# -eq 0 ]]; then usage; exit 1; fi

while [[ $# -gt 0 ]]; do
  case "$1" in
    -g|--genome-dir) FASTA_DIR="$2"; shift 2;;
    -r|--rnaseq-dir) RNA_SEQ_DIR="$2"; shift 2;;
    -o|--outdir)     OUTDIR="$2"; shift 2;;

    -i|--image)      SINGULARITY_IMAGE="$2"; shift 2;;
    -t|--threads)    THREADS="$2"; shift 2;;

    --suffix)        GENOME_SUFFIX="$2"; shift 2;;

    --no-rna-bind)   BIND_RNA=false; shift;;
    --cleanup)       NO_CLEANUP=false; shift;;

    --no-fungus)     FUNGUS=false; shift;;
    --no-ab-initio)  AB_INITIO=false; shift;;
    --no-gff3)       GFF3=false; shift;;

    -h|--help) usage; exit 0;;
    *) die "Unknown option: $1";;
  esac
done
[[ -z "${THREADS}" ]] || positive_int THREADS "${THREADS}"

# -------------------------
# validate
# -------------------------
[[ -d "$FASTA_DIR" ]]    || die "Genome dir not found"
[[ -d "$RNA_SEQ_DIR" ]]  || die "RNA-seq dir not found"
[[ -n "$OUTDIR" ]]       || die "Missing --outdir"
[[ -f "$SINGULARITY_IMAGE" ]] || die "Singularity image not found"

command -v singularity >/dev/null || die "singularity not found"

mkdir -p "$OUTDIR"

FASTA_DIR=$(realpath "$FASTA_DIR")
RNA_SEQ_DIR=$(realpath "$RNA_SEQ_DIR")
OUTDIR=$(realpath "$OUTDIR")
SINGULARITY_IMAGE=$(realpath "$SINGULARITY_IMAGE")

echo "=============================="
echo "BRaKER3 batch configuration"
echo "  genome_dir   : $FASTA_DIR"
echo "  rnaseq_dir   : $RNA_SEQ_DIR"
echo "  outdir       : $OUTDIR"
echo "  image        : $SINGULARITY_IMAGE"
echo "  threads      : $THREADS"
echo "  suffix       : $GENOME_SUFFIX"
echo "=============================="

# -------------------------
# main loop (top-level only)
# -------------------------
shopt -s nullglob
genomes=( "$FASTA_DIR"/*"$GENOME_SUFFIX" )
(( ${#genomes[@]} > 0 )) || die "No genome FASTA found"

for GENOME_FILE in "${genomes[@]}"; do
  SAMPLE=$(basename "$GENOME_FILE" "$GENOME_SUFFIX")
  WORKDIR="$OUTDIR/$SAMPLE"
  mkdir -p "$WORKDIR"

  echo "=== Running BRaKER3 for $SAMPLE ==="

  BIND_OPTS=(-B "$FASTA_DIR:$FASTA_DIR" -B "$OUTDIR:$OUTDIR")
  CONTAINER_RNA_DIR="$RNA_SEQ_DIR"
  if [[ "$BIND_RNA" == true ]]; then
    BIND_OPTS+=(-B "$RNA_SEQ_DIR:/RNA_seq")
    CONTAINER_RNA_DIR="/RNA_seq"
  fi

  BRAKER_OPTS=()
  $FUNGUS && BRAKER_OPTS+=(--fungus)
  $AB_INITIO && BRAKER_OPTS+=(--AUGUSTUS_ab_initio)
  $GFF3 && BRAKER_OPTS+=(--gff3)
  $NO_CLEANUP && BRAKER_OPTS+=(--nocleanup)

  singularity exec \
    "${BIND_OPTS[@]}" \
    "$SINGULARITY_IMAGE" \
    braker.pl \
      --genome="$GENOME_FILE" \
      --rnaseq_sets_ids="$SAMPLE" \
      --rnaseq_sets_dir="$CONTAINER_RNA_DIR" \
      --species="$SAMPLE" \
      --threads "$THREADS" \
      --workingdir="$WORKDIR" \
      "${BRAKER_OPTS[@]}"

  echo "=== Finished $SAMPLE ==="
done
