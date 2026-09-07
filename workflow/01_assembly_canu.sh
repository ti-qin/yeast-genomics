#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/lib/common.sh"
validate_value_options "--canu --genome-size --outdir --reads -g -o -r" "$@"

#######################################
# defaults
#######################################
READSDIR=""
OUTDIR="./result/assembly/canu"
CANU_CMD="canu"
GENOME_SIZE="12m"
OVERWRITE=false

#######################################
# help
#######################################
usage() {
cat <<EOF
Usage: $(basename "$0") -r READSDIR [options]

Required:
  -r, --reads DIR         Directory containing ONT fastq files
                          (e.g. G12357.fastq.gz)

Optional:
  -o, --outdir DIR        Output directory
                          (default: ./result/assembly/canu)
  -g, --genome-size STR   Genome size for Canu
                          (default: 12m)
      --canu PATH         Canu executable (default: canu in PATH)
      --overwrite         Overwrite existing sample directories
  -h, --help              Show this help message

Behavior:
  - Samples are processed sequentially (no parallel jobs)
  - Each fastq file is treated as one sample
  - Output layout:
      <outdir>/<sample>/
        ├── <sample>.contigs.fasta
        ├── <sample>.canu.log
        └── ...

Example:
  $(basename "$0") -r ./rawdata/ont -o ./result/assembly/canu -g 12m

nohup example:
  nohup $(basename "$0") \\
    -r ./rawdata/ont \\
    -o ./result/assembly/canu \\
    -g 12m \\
    > canu.batch.log 2>&1 &
EOF
}

#######################################
# parse arguments
#######################################
while [[ $# -gt 0 ]]; do
  case "$1" in
    -r|--reads)
      READSDIR="$2"; shift 2;;
    -o|--outdir)
      OUTDIR="$2"; shift 2;;
    -g|--genome-size)
      GENOME_SIZE="$2"; shift 2;;
    --canu)
      CANU_CMD="$2"; shift 2;;
    --overwrite)
      OVERWRITE=true; shift;;
    -h|--help)
      usage; exit 0;;
    *)
      echo "ERROR: Unknown option: $1"
      usage; exit 1;;
  esac
done

#######################################
# checks
#######################################
if [[ -z "$READSDIR" ]]; then
  echo "ERROR: --reads is required"
  usage
  exit 1
fi

command -v "$CANU_CMD" >/dev/null 2>&1 || {
  echo "ERROR: canu not found: $CANU_CMD"
  exit 1
}

mkdir -p "$OUTDIR"

#######################################
# find fastq files
#######################################
shopt -s nullglob
FASTQ_FILES=("$READSDIR"/*.fastq.gz "$READSDIR"/*.fq.gz "$READSDIR"/*.fastq "$READSDIR"/*.fq)
shopt -u nullglob

if [[ ${#FASTQ_FILES[@]} -eq 0 ]]; then
  echo "ERROR: No fastq files found in $READSDIR"
  exit 1
fi

#######################################
# summary
#######################################
echo "=== Canu serial batch assembly started: $(date) ==="
echo "Reads dir:    $READSDIR"
echo "Output dir:   $OUTDIR"
echo "Genome size:  $GENOME_SIZE"
echo "Canu cmd:     $CANU_CMD"
echo "Overwrite:    $OVERWRITE"
echo "Samples:      ${#FASTQ_FILES[@]}"
echo "==================================================="

#######################################
# main loop (serial)
#######################################
for fq in "${FASTQ_FILES[@]}"; do
  fname=$(basename "$fq")
  sample=${fname%%.fastq.gz}
  sample=${sample%%.fq.gz}
  sample=${sample%%.fastq}
  sample=${sample%%.fq}

  sample_outdir="$OUTDIR/$sample"
  sample_log="$sample_outdir/$sample.canu.log"

  if [[ -d "$sample_outdir" && "$OVERWRITE" != true ]]; then
    echo "⏭️  Skip $sample (directory exists)"
    continue
  fi

  rm -rf "$sample_outdir"
  mkdir -p "$sample_outdir"

  echo
  echo "🧬 [$(date '+%F %T')] START $sample"
  echo "    Reads : $fq"
  echo "    Outdir: $sample_outdir"

  "$CANU_CMD" \
    -p "$sample" \
    -d "$sample_outdir" \
    genomeSize="$GENOME_SIZE" \
    -nanopore-raw "$fq" \
    > "$sample_log" 2>&1

  echo "✅ [$(date '+%F %T')] DONE $sample"
done

echo
echo "🎉 All Canu jobs finished: $(date)"
