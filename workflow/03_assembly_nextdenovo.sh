#!/bin/bash
# ===============================================================
# Batch long-read assembly using nextDenovo (parameter mode)
# + optional clean: remove contigs with type:s:loop from fasta
# ===============================================================

set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/lib/common.sh"
validate_value_options "--genome --jobs --reads --threads --workdir -g -j -r -t -w" "$@"

READS_DIR="./reads"
WORKDIR="./"
THREADS=16
PARALLEL_JOBS=2
GENOME_SIZE=12000000

NEXTDENOVO_CMD="nextDenovo"

# new
CLEAN=false

usage() {
cat <<EOF
Usage: $(basename "$0") [options]

Options:
  -r, --reads DIR        reads directory (default: ./reads)
  -w, --workdir DIR      work/output directory (default: ./)
  -t, --threads INT      threads for minimap2 (default: 16)
  -j, --jobs INT         concurrent samples AND internal jobs per sample (default: 2)
  -g, --genome INT       genome size bp (default: 12000000)
      --clean            after assembly, remove contigs whose header matches:
                         'type:s:loop'
                         output: <workdir>/<sample>.nextdenovo.clean.fasta
  -h, --help             show help

Example:
  $(basename "$0") -r ./raw_data/ont -w ./result/assembly/denovo -t 16 -j 2 -g 12000000 --clean

nohup:
  nohup $(basename "$0") -r ./raw_data/ont -w ./result/assembly/denovo -t 16 -j 2 -g 12000000 --clean \
    > nextdenovo.nohup.log 2>&1 &
EOF
}

# ---------------------------
# parse args
# ---------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    -r|--reads)   READS_DIR="$2"; shift 2 ;;
    -w|--workdir) WORKDIR="$2"; shift 2 ;;
    -t|--threads) THREADS="$2"; shift 2 ;;
    -j|--jobs)    PARALLEL_JOBS="$2"; shift 2 ;;
    -g|--genome)  GENOME_SIZE="$2"; shift 2 ;;
    --clean)      CLEAN=true; shift ;;
    -h|--help)    usage; exit 0 ;;
    *) echo "ERROR: Unknown option: $1"; usage; exit 1 ;;
  esac
done
[[ -z "${THREADS}" ]] || positive_int THREADS "${THREADS}"
[[ -z "${PARALLEL_JOBS}" ]] || positive_int PARALLEL_JOBS "${PARALLEL_JOBS}"

# ---------------------------
# basic checks
# ---------------------------
if ! command -v "$NEXTDENOVO_CMD" >/dev/null 2>&1; then
  echo "❌ ERROR: nextDenovo not found in PATH"
  echo "   Please activate the environment containing nextDenovo"
  exit 1
fi

if [[ "$CLEAN" == true ]]; then
  if ! command -v seqkit >/dev/null 2>&1; then
    echo "❌ ERROR: seqkit not found but --clean requires seqkit"
    exit 1
  fi
fi

mkdir -p "$WORKDIR"
WORKDIR=$(realpath "$WORKDIR")
READS_DIR=$(realpath "$READS_DIR")
PIDS=()

echo "=== nextDenovo batch started: $(date) ==="
echo "READS_DIR:     $READS_DIR"
echo "WORKDIR:       $WORKDIR"
echo "THREADS:       $THREADS"
echo "PARALLEL_JOBS: $PARALLEL_JOBS"
echo "GENOME_SIZE:   $GENOME_SIZE"
echo "nextDenovo:    $(command -v nextDenovo)"
echo "CLEAN:         $CLEAN"
echo "========================================="

# ---------------------------
# find reads
# ---------------------------
shopt -s nullglob
READ_FILES=("$READS_DIR"/*.fastq.gz)
shopt -u nullglob

if [[ ${#READ_FILES[@]} -eq 0 ]]; then
  echo "❌ No fastq.gz files found in $READS_DIR"
  exit 1
fi

# ---------------------------
# main loop
# ---------------------------
for READ_FILE in "${READ_FILES[@]}"; do
  PREFIX=$(basename "$READ_FILE" .fastq.gz)
  SAMPLE_DIR="$WORKDIR/$PREFIX"
  mkdir -p "$SAMPLE_DIR"

  echo "🚀 Starting nextDenovo for $PREFIX"
  echo "   Reads:       $READ_FILE"
  echo "   Workdir:     $SAMPLE_DIR"
  echo "   Genome size: $GENOME_SIZE bp"

  # input.fofn
  FOFN="$SAMPLE_DIR/input.fofn"
  FOFN_ABS=$(realpath "$FOFN")
  READ_ABS=$(realpath "$READ_FILE")
  echo "$READ_ABS" > "$FOFN"

  # run.cfg
  CFG=$(realpath "$SAMPLE_DIR/run.cfg")
  cat > "$CFG" <<EOF
[General]
job_type = local
job_prefix = nextDenovo
task = all
rewrite = yes
deltmp = yes
rerun = 3
parallel_jobs = $PARALLEL_JOBS
input_type = raw
read_type = ont
input_fofn = $FOFN_ABS
workdir = $SAMPLE_DIR

[correct_option]
read_cutoff = 1k
genome_size = $GENOME_SIZE
pa_correction = 2
sort_options = -m 1g -t 2
minimap2_options_raw = -t $THREADS
correction_options = -p 15

[assemble_option]
minimap2_options_cns = -t $THREADS
nextgraph_options = -a 1
EOF

  (
    echo "📊 Running nextDenovo for $PREFIX..."
    "$NEXTDENOVO_CMD" "$CFG"

    FINAL_FASTA="$SAMPLE_DIR/03.ctg_graph/nd.asm.fasta"
    OUT_FASTA="$WORKDIR/${PREFIX}.nextdenovo.raw.fasta"
    OUT_CLEAN="$WORKDIR/${PREFIX}.nextdenovo.fasta"

    if [[ -s "$FINAL_FASTA" ]]; then
      cp -f "$FINAL_FASTA" "$OUT_FASTA"
      echo "✅ Finished: $OUT_FASTA"

      if [[ "$CLEAN" == true ]]; then
        echo "🧹 Cleaning loop contigs: $OUT_CLEAN"
        # remove contigs whose header contains type:s:loop
        seqkit grep -r -n -v -p 'type:s:loop' "$OUT_FASTA" > "$OUT_CLEAN"
        rm -f "$OUT_FASTA"
        echo "✅ Cleaned: $OUT_CLEAN"
      fi
    else
      echo "❌ No assembly fasta found for $PREFIX"
      echo "   Expected: $FINAL_FASTA"
      exit 1
    fi
  ) &
  PIDS+=("$!")


  while (( $(jobs -r | wc -l) >= PARALLEL_JOBS )); do
    echo "⏳ Waiting for slots... ($(jobs -r | wc -l) jobs running)"
    sleep 3
  done
done

FAILED=0
for pid in "${PIDS[@]}"; do
  wait "$pid" || FAILED=1
done
(( FAILED == 0 )) || { echo "ERROR: One or more samples failed" >&2; exit 1; }
echo "🎉 All nextDenovo assembly jobs finished!"
