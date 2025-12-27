#!/usr/bin/env bash
set -euo pipefail

#######################################
# defaults
#######################################
READSDIR=""
OUTDIR="."
FLYE_CMD="flye"
THREADS_PER_JOB=32
MAX_JOBS=""
FORCE_OVERWRITE=false

#######################################
# help
#######################################
usage() {
cat << EOF
Usage: $(basename "$0") -r READSDIR [options]

Required:
  -r, --reads DIR         Directory containing fastq files

Optional:
  -o, --outdir DIR        Output directory (default: .)
  -t, --threads INT       Threads per Flye job (default: 32)
  -j, --jobs INT          Max concurrent jobs (default: auto)
      --flye PATH         Flye executable (default: flye in PATH)
      --overwrite         Overwrite existing output directories
  -h, --help              Show this help message

Example:
  $(basename "$0") -r ../reads -o flye_out -t 32 -j 4 --overwrite
nohup example:
  nohup $(basename "$0") -r ./raw_data/ont -o ./result/assembly/flye -t 32 -j 4 --overwrite > flye.nohup.log 2>&1 &
EOF
}

#######################################
# parse arguments
#######################################
while [[ $# -gt 0 ]]; do
  case "$1" in
    -r|--reads)
      READSDIR="$2"
      shift 2
      ;;
    -o|--outdir)
      OUTDIR="$2"
      shift 2
      ;;
    -t|--threads)
      THREADS_PER_JOB="$2"
      shift 2
      ;;
    -j|--jobs)
      MAX_JOBS="$2"
      shift 2
      ;;
    --flye)
      FLYE_CMD="$2"
      shift 2
      ;;
    --overwrite)
      FORCE_OVERWRITE=true
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "ERROR: Unknown option: $1"
      usage
      exit 1
      ;;
  esac
done

#######################################
# sanity checks
#######################################
if [[ -z "$READSDIR" ]]; then
  echo "ERROR: --reads is required"
  usage
  exit 1
fi

mkdir -p "$OUTDIR"

#######################################
# detect fastq files
#######################################
shopt -s nullglob
FASTQ_FILES=("$READSDIR"/*.fastq.gz "$READSDIR"/*.fq.gz "$READSDIR"/*.fastq "$READSDIR"/*.fq)
shopt -u nullglob

if [[ ${#FASTQ_FILES[@]} -eq 0 ]]; then
  echo "ERROR: No fastq files found in $READSDIR"
  exit 1
fi

#######################################
# auto-calc MAX_JOBS if not provided
#######################################
CORES=$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 1)

if [[ -z "$MAX_JOBS" ]]; then
  MAX_JOBS=$(( CORES / THREADS_PER_JOB ))
  (( MAX_JOBS < 1 )) && MAX_JOBS=1
fi

#######################################
# summary
#######################################
echo "=== Flye batch started: $(date) ==="
echo "Reads dir:     $READSDIR"
echo "Output dir:    $OUTDIR"
echo "Flye cmd:      $FLYE_CMD"
echo "Threads/job:   $THREADS_PER_JOB"
echo "Max jobs:      $MAX_JOBS"
echo "Overwrite:     $FORCE_OVERWRITE"
echo "Detected CPUs: $CORES"
echo "==================================="

#######################################
# helpers
#######################################
running_jobs() { jobs -rp | wc -l; }

wait_for_slot() {
  while (( $(running_jobs) >= MAX_JOBS )); do
    sleep 3
  done
}

FAILURES_FILE="$(mktemp --tmpdir flye_failures.XXXXXX)"
trap 'echo "Interrupted, waiting jobs..."; wait; exit 1' INT TERM

#######################################
# main loop
#######################################
for fq in "${FASTQ_FILES[@]}"; do
  fname=$(basename "$fq")
  sample=${fname%%.fastq.gz}
  sample=${sample%%.fq.gz}
  sample=${sample%%.fastq}
  sample=${sample%%.fq}

  sample_outdir="$OUTDIR/$sample"
  sample_log="$sample_outdir/$sample.log"

  if [[ -d "$sample_outdir" && "$FORCE_OVERWRITE" != true ]]; then
    echo "Skip $sample (exists)"
    continue
  fi

  rm -rf "$sample_outdir"
  mkdir -p "$sample_outdir"

  CMD=(
    "$FLYE_CMD"
    --nano-raw "$fq"
    --out-dir "$sample_outdir"
    --threads "$THREADS_PER_JOB"
    --iterations 3
  )
  
  wait_for_slot

  (
    echo "[$(date)] START $sample"
    echo "CMD: ${CMD[*]}"
    if "${CMD[@]}" > "$sample_log" 2>&1; then
      echo "[$(date)] DONE $sample" >> "$sample_log"
      FINAL_FASTA="${sample_outdir}/assembly.fasta"
      FINAL_INFO="${sample_outdir}/assembly_info.txt"
      OUT_FASTA="${OUTDIR}/${sample}.flye.fasta"
      OUT_INFO="${OUTDIR}/${sample}.flye_info.txt"
      if [[ ! -s "$FINAL_FASTA" || ! -s "$FINAL_INFO" ]]; then
        echo "[$(date)] ERROR: Flye output missing for $sample" >> "$sample_log"
        exit 1
      fi
      cp -f "$FINAL_FASTA" "$OUT_FASTA"
      cp -f "$FINAL_INFO" "$OUT_INFO"
      echo "[$(date)] Exported Flye results" >> "$sample_log"
    else
      echo "[$(date)] FAIL $sample" >> "$sample_log"
      echo "$sample" >> "$FAILURES_FILE"
    fi
  ) &
done

wait

#######################################
# report
#######################################
if [[ -s "$FAILURES_FILE" ]]; then
  echo "Some samples failed:"
  cat "$FAILURES_FILE"
  exit 2
else
  echo "🎉 All samples completed successfully"
fi
