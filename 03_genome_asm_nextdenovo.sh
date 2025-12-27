#!/bin/bash
# ===============================================================
# Batch long-read assembly using nextDenovo (parameter mode)
# ===============================================================

set -euo pipefail

READS_DIR="./reads"
WORKDIR="./"
THREADS=16
PARALLEL_JOBS=2
GENOME_SIZE=12000000

NEXTDENOVO_CMD="nextDenovo"  

usage() {
cat <<EOF
Usage: $(basename "$0") [options]

Options:
  -r, --reads DIR        reads directory (default: ./raw_data/ont)
  -w, --workdir DIR      work/output directory (default: ./result/assembly/denovo)
  -t, --threads INT      threads for minimap2 (default: 16)
  -j, --jobs INT         parallel samples (default: 2)
  -g, --genome INT       genome size bp (default: 12000000)
  -h, --help             show help

Example:
  $(basename "$0") -r ./raw_data/ont -w ./result/assembly/denovo -t 16 -j 2 -g 12000000

nohup:
  nohup $(basename "$0") -r ./raw_data/ont -w ./result/assembly/denovo -t 16 -j 2 -g 12000000 \
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
    -h|--help)    usage; exit 0 ;;
    *) echo "ERROR: Unknown option: $1"; usage; exit 1 ;;
  esac
done

# ---------------------------
# basic checks
# ---------------------------
if ! command -v "$NEXTDENOVO_CMD" >/dev/null 2>&1; then
  echo "❌ ERROR: nextDenovo not found in PATH"
  echo "   Please activate the environment containing nextDenovo"
  exit 1
fi

mkdir -p "$WORKDIR"

echo "=== nextDenovo batch started: $(date) ==="
echo "READS_DIR:     $READS_DIR"
echo "WORKDIR:       $WORKDIR"
echo "THREADS:       $THREADS"
echo "PARALLEL_JOBS: $PARALLEL_JOBS"
echo "GENOME_SIZE:   $GENOME_SIZE"
echo "nextDenovo:    $(command -v nextDenovo)"
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

# 保留你的原逻辑：复制 reads 到样本目录
#   cp -f "$READ_FILE" "$SAMPLE_DIR/"
#   LOCAL_READ_FILE=$(basename "$READ_FILE")

  echo "🚀 Starting nextDenovo for $PREFIX"
  echo "   Reads:       $READ_FILE"
  echo "   Workdir:     $SAMPLE_DIR"
  echo "   Genome size: $GENOME_SIZE bp"

  # input.fofn（写本地文件名）
  FOFN="$SAMPLE_DIR/input.fofn"
  FOFN_ABS=$(realpath "$FOFN")
  READ_ABS=$(realpath "$READ_FILE")
  echo "$READ_ABS" > "$FOFN"

  # run.cfg（结构与你原脚本一致）
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
workdir = .

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
    nextDenovo "$CFG"

    FINAL_FASTA="$SAMPLE_DIR/03.ctg_graph/nd.asm.fasta"
    if [[ -f "$FINAL_FASTA" ]]; then
      cp -f "$FINAL_FASTA" "$WORKDIR/${PREFIX}.nextdenovo.fasta"
      echo "✅ Finished: $WORKDIR/${PREFIX}.nextdenovo.fasta"
    else
      echo "❌ No assembly fasta found for $PREFIX"
      echo "   Expected: $FINAL_FASTA"
    fi
  ) &

  # 控制样本并发数（保持你原来的方式）
  while (( $(jobs -r | wc -l) >= PARALLEL_JOBS )); do
    echo "⏳ Waiting for slots... ($(jobs -r | wc -l) jobs running)"
    sleep 3
  done
done

wait
echo "🎉 All nextDenovo assembly jobs finished!"
