#!/usr/bin/env bash
set -euo pipefail

#######################################
# defaults (can be overridden by args)
#######################################
READSDIR=""
OUTDIR="."
NECAT_CMD="necat.pl"
MAX_JOBS=4

# NECAT options
GENOME_SIZE=12000000
THREADS=32
MIN_READ_LENGTH=1000
PREP_OUTPUT_COVERAGE=40
CNS_OUTPUT_COVERAGE=30
NUM_ITER=2
POLISH_CONTIGS=true

# NECAT option strings (allow override)
OVLP_FAST_OPTIONS='-n 500 -z 20 -b 2000 -e 0.5 -j 0 -u 1 -a 1000'
OVLP_SENSITIVE_OPTIONS='-n 500 -z 10 -e 0.5 -j 0 -u 1 -a 1000'
CNS_FAST_OPTIONS='-a 2000 -x 4 -y 12 -l 1000 -e 0.5 -p 0.8 -u 0'
CNS_SENSITIVE_OPTIONS='-a 2000 -x 4 -y 12 -l 1000 -e 0.5 -p 0.8 -u 0'
TRIM_OVLP_OPTIONS='-n 100 -z 10 -b 2000 -e 0.5 -j 1 -u 1 -a 400'
ASM_OVLP_OPTIONS='-n 100 -z 10 -b 2000 -e 0.5 -j 1 -u 0 -a 400'

# behavior
FORCE_OVERWRITE=false
KEEP_CONFIG=false   # default clean temp config/readlist, set true to keep

usage() {
cat << EOF
Usage: $(basename "$0") -r READSDIR [options]

Required:
  -r, --reads DIR            Directory containing fastq(.gz)/fq(.gz)

Optional:
  -o, --outdir DIR           Output directory (default: .)
  -n, --necat PATH           necat.pl path or command (default: necat.pl)
  -j, --jobs INT             Max concurrent samples (default: 4)
  -t, --threads INT          Threads per sample (default: 32)
  -g, --genome-size INT      Genome size (default: 12000000)
      --min-read-len INT     MIN_READ_LENGTH (default: 1000)
      --prep-cov INT         PREP_OUTPUT_COVERAGE (default: 40)
      --cns-cov INT          CNS_OUTPUT_COVERAGE (default: 30)
      --num-iter INT         NUM_ITER (default: 2)
      --polish true|false    POLISH_CONTIGS (default: true)

      --ovlp-fast "..."      OVLP_FAST_OPTIONS (default built-in)
      --ovlp-sensitive "..." OVLP_SENSITIVE_OPTIONS
      --cns-fast "..."       CNS_FAST_OPTIONS
      --cns-sensitive "..."  CNS_SENSITIVE_OPTIONS
      --trim-ovlp "..."      TRIM_OVLP_OPTIONS
      --asm-ovlp "..."       ASM_OVLP_OPTIONS

      --overwrite            Overwrite existing sample output dirs
      --keep-config          Keep generated readlist/config files
  -h, --help                 Show help

Example:
  $(basename "$0") -r ../reads -o ./necat_out -g 12000000 -t 32 -j 4 --overwrite

nohup example:
  nohup $(basename "$0") -r ./raw_data/ont -o ./result/assembly/necat -t 32 -j 4 --overwrite > necat.nohup.log 2>&1 &
EOF
}

#######################################
# parse args
#######################################
if [[ $# -eq 0 ]]; then
  usage
  exit 1
fi

while [[ $# -gt 0 ]]; do
  case "$1" in
    -r|--reads) READSDIR="$2"; shift 2;;
    -o|--outdir) OUTDIR="$2"; shift 2;;
    -n|--necat) NECAT_CMD="$2"; shift 2;;
    -j|--jobs) MAX_JOBS="$2"; shift 2;;
    -t|--threads) THREADS="$2"; shift 2;;
    -g|--genome-size) GENOME_SIZE="$2"; shift 2;;

    --min-read-len) MIN_READ_LENGTH="$2"; shift 2;;
    --prep-cov) PREP_OUTPUT_COVERAGE="$2"; shift 2;;
    --cns-cov) CNS_OUTPUT_COVERAGE="$2"; shift 2;;
    --num-iter) NUM_ITER="$2"; shift 2;;
    --polish) POLISH_CONTIGS="$2"; shift 2;;

    --ovlp-fast) OVLP_FAST_OPTIONS="$2"; shift 2;;
    --ovlp-sensitive) OVLP_SENSITIVE_OPTIONS="$2"; shift 2;;
    --cns-fast) CNS_FAST_OPTIONS="$2"; shift 2;;
    --cns-sensitive) CNS_SENSITIVE_OPTIONS="$2"; shift 2;;
    --trim-ovlp) TRIM_OVLP_OPTIONS="$2"; shift 2;;
    --asm-ovlp) ASM_OVLP_OPTIONS="$2"; shift 2;;

    --overwrite) FORCE_OVERWRITE=true; shift;;
    --keep-config) KEEP_CONFIG=true; shift;;

    -h|--help) usage; exit 0;;
    *) echo "ERROR: Unknown option: $1"; usage; exit 1;;
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

# find fastq files
shopt -s nullglob
FASTQ_FILES=("$READSDIR"/*.fastq.gz "$READSDIR"/*.fq.gz "$READSDIR"/*.fastq "$READSDIR"/*.fq)
shopt -u nullglob

if [[ ${#FASTQ_FILES[@]} -eq 0 ]]; then
  echo "❌ ERROR: can not find fastq/fq in: $READSDIR"
  exit 1
fi

echo "=== NECAT parallel batch assembly started: $(date) ==="
echo "Reads dir:   $READSDIR"
echo "Outdir:      $OUTDIR"
echo "NECAT:       $NECAT_CMD"
echo "Samples:     ${#FASTQ_FILES[@]}"
echo "Threads/job: $THREADS"
echo "Max jobs:    $MAX_JOBS"
echo "Overwrite:   $FORCE_OVERWRITE"
echo "Keep config: $KEEP_CONFIG"
echo "Genome size: $GENOME_SIZE"
echo "=============================================="

#######################################
# concurrency helpers
#######################################
wait_for_jobs() {
  while (( $(jobs -rp | wc -l) >= MAX_JOBS )); do
    sleep 5
  done
}

FAILURES_FILE="$(mktemp --tmpdir necat_failures.XXXXXX)"
trap 'echo "Interrupted, waiting jobs..."; wait; exit 1' INT TERM

#######################################
# main loop
#######################################
for READS in "${FASTQ_FILES[@]}"; do
  BASENAME=$(basename "$READS")
  SAMPLE=${BASENAME%%.fastq.gz}
  SAMPLE=${SAMPLE%%.fq.gz}
  SAMPLE=${SAMPLE%%.fastq}
  SAMPLE=${SAMPLE%%.fq}

  PREFIX="${SAMPLE}"
  SAMPLE_DIR="${OUTDIR}/${SAMPLE}"
  mkdir -p "$SAMPLE_DIR"

  # overwrite behavior
  if [[ -e "$SAMPLE_DIR/${PREFIX}" || -e "$SAMPLE_DIR/${PREFIX}.config.txt" || -d "$SAMPLE_DIR/${PREFIX}" ]]; then
    if [[ "$FORCE_OVERWRITE" != true ]]; then
      echo "Skip $SAMPLE (existing output in $SAMPLE_DIR). Use --overwrite to rerun."
      continue
    fi
    # cautious cleanup (remove only this sample dir contents)
    rm -rf "$SAMPLE_DIR"
    mkdir -p "$SAMPLE_DIR"
  fi

  SAMPLE_LOG="${SAMPLE_DIR}/${SAMPLE}.necat.log"
  READLIST="${SAMPLE_DIR}/reads_${PREFIX}.txt"
  CONFIG="${SAMPLE_DIR}/${PREFIX}.config.txt"

  wait_for_jobs

  (
    echo "🧬 [$(date '+%F %T')] START: ${SAMPLE}"
    echo "Reads: $READS"
    echo "Log:   $SAMPLE_LOG"

    echo "$READS" > "$READLIST"

    cat > "$CONFIG" <<EOF
PROJECT=${SAMPLE_DIR}
ONT_READ_LIST=${READLIST}
GENOME_SIZE=${GENOME_SIZE}
THREADS=${THREADS}
MIN_READ_LENGTH=${MIN_READ_LENGTH}
PREP_OUTPUT_COVERAGE=${PREP_OUTPUT_COVERAGE}
OVLP_FAST_OPTIONS=${OVLP_FAST_OPTIONS}
OVLP_SENSITIVE_OPTIONS=${OVLP_SENSITIVE_OPTIONS}
CNS_FAST_OPTIONS=${CNS_FAST_OPTIONS}
CNS_SENSITIVE_OPTIONS=${CNS_SENSITIVE_OPTIONS}
TRIM_OVLP_OPTIONS=${TRIM_OVLP_OPTIONS}
ASM_OVLP_OPTIONS=${ASM_OVLP_OPTIONS}
NUM_ITER=${NUM_ITER}
CNS_OUTPUT_COVERAGE=${CNS_OUTPUT_COVERAGE}
CLEANUP=1
USE_GRID=false
SMALL_MEMORY=0
POLISH_CONTIGS=${POLISH_CONTIGS}
EOF

    # run NECAT (redirect ALL output to sample log)
    {
      echo "👉 [${SAMPLE}] correct"
      "${NECAT_CMD}" correct "$CONFIG"

      echo "👉 [${SAMPLE}] assemble"
      "${NECAT_CMD}" assemble "$CONFIG"

      echo "👉 [${SAMPLE}] bridge"
      "${NECAT_CMD}" bridge "$CONFIG"

       echo "✅ [$(date '+%F %T')] DONE: ${SAMPLE}"
    }> "$SAMPLE_LOG" 2>&1

    # cleanup configs unless requested
    if [[ "$KEEP_CONFIG" != true ]]; then
      rm -f "$READLIST" "$CONFIG"
    fi

    FINAL_FASTA="$SAMPLE_DIR/6-bridge_contigs/polished_contigs.fasta"
    if [[ -f "$FINAL_FASTA" ]]; then
      cp -f "$FINAL_FASTA" "$OUTDIR/${PREFIX}.necat.fasta"
      echo "✅ Finished: $OUTDIR/${PREFIX}.necat.fasta"
    else
      echo "❌ No assembly fasta found for $PREFIX"
      echo "   Expected: $FINAL_FASTA"
    fi

  ) || {
    echo "❌ [$(date '+%F %T')] FAIL: ${SAMPLE}"
    echo "$SAMPLE" >> "$FAILURES_FILE"
  } &
done

wait

#######################################
# report
#######################################
if [[ -s "$FAILURES_FILE" ]]; then
  echo "=== Some samples failed ==="
  while read -r s; do
    echo "FAILED: $s  (see ${OUTDIR}/${s}/${s}.necat.log)"
  done < "$FAILURES_FILE"
  exit 2
else
  echo "🎉 All done: $(date)"
fi
