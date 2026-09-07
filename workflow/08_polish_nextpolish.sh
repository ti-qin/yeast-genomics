#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/lib/common.sh"
validate_value_options "--assembly-dir --jobs --nextpolish --outdir --pattern --r1-regex --r2-regex --reads-dir --suffix --threads -a -j -o -r -t" "$@"

# ===============================================================
# Batch short-read polishing using NextPolish (CLI tool)
# ===============================================================

# -----------------------------
# defaults
# -----------------------------
ASSEMBLY_DIR="./racon_res"
READS_DIR="../../illumina"
WORKDIR="./nextpo_res"
THREADS=32
PARALLEL_JOBS=2

# how to find assemblies
ASM_GLOB="*.recon3.fasta"
ASM_SUFFIX=".recon3.fasta"

# how to find reads (regex suffix)
R1_REGEX='[_.-]1\.fq\.gz$'
R2_REGEX='[_.-]2\.fq\.gz$'

# nextpolish executable
NEXTPOLISH="nextPolish"

# behavior
COPY_READS=false     # default: use symlink (saves space)
KEEP_SAMPLE_DIR=true # default: keep per-sample folders

# -----------------------------
# help
# -----------------------------
usage() {
  cat <<'EOF'
Usage:
  08_polish_nextpolish.sh -a ASSEMBLY_DIR -r READS_DIR -o OUTDIR [options]

Required/Typical:
  -a, --assembly-dir DIR     Directory containing assemblies (default: ./racon_res)
  -r, --reads-dir DIR        Directory containing Illumina reads (default: ../../illumina)
  -o, --outdir DIR           Output directory (default: ./nextpo_res)

Options:
  -t, --threads INT          Threads per NextPolish job (default: 32)
  -j, --jobs INT             Max concurrent samples (default: 2)

  --pattern GLOB             Assembly filename glob under ASSEMBLY_DIR
                             (default: *.recon3.fasta)
  --suffix STR               Suffix to strip for PREFIX (default: .recon3.fasta)

  --r1-regex REGEX           Regex to match R1 fastq.gz (default: [_.-]1\.fq\.gz$)
  --r2-regex REGEX           Regex to match R2 fastq.gz (default: [_.-]2\.fq\.gz$)

  --nextpolish PATH          Path to nextPolish executable
                             (default: nextPolish in PATH)

  --copy-reads               Copy reads into per-sample folder (default: symlink)
  --no-keep-sample-dir       Remove per-sample workdir after success (keeps final fasta)

  -h, --help                 Show this help

Naming rules:
  PREFIX is inferred from assembly filename by stripping --suffix.
  Reads are searched within READS_DIR by:
    .*/${PREFIX}.*<R1_REGEX> and .*/${PREFIX}.*<R2_REGEX>
  (exactly one match per mate is required)

Examples:
  # default pattern/suffix
  ./08_polish_nextpolish.sh -a ./racon_res -r ../../illumina -o ./nextpo_res -t 32 -j 2

  # different assembly suffix/pattern
  workflow/08_polish_nextpolish -a ./asm -r ./reads -o ./np \
    --pattern "*.fasta" --suffix ".fasta" -t 48 -j 4

  # specify nextPolish binary
  workflow/08_polish_nextpolish -a ./asm -r ./reads -o ./np --nextpolish "$(which nextPolish)"

EOF
}

die() { echo "❌ $*" >&2; exit 1; }

# -----------------------------
# arg parsing
# -----------------------------
if [[ $# -eq 0 ]]; then usage; exit 1; fi

while [[ $# -gt 0 ]]; do
  case "$1" in
    -a|--assembly-dir) ASSEMBLY_DIR="$2"; shift 2;;
    -r|--reads-dir)    READS_DIR="$2"; shift 2;;
    -o|--outdir)       WORKDIR="$2"; shift 2;;
    -t|--threads)      THREADS="$2"; shift 2;;
    -j|--jobs)         PARALLEL_JOBS="$2"; shift 2;;

    --pattern)         ASM_GLOB="$2"; shift 2;;
    --suffix)          ASM_SUFFIX="$2"; shift 2;;

    --r1-regex)        R1_REGEX="$2"; shift 2;;
    --r2-regex)        R2_REGEX="$2"; shift 2;;

    --nextpolish)      NEXTPOLISH="$2"; shift 2;;

    --copy-reads)      COPY_READS=true; shift 1;;
    --no-keep-sample-dir) KEEP_SAMPLE_DIR=false; shift 1;;

    -h|--help) usage; exit 0;;
    *) die "Unknown option: $1 (use -h for help)";;
  esac
done
[[ -z "${THREADS}" ]] || positive_int THREADS "${THREADS}"
[[ -z "${PARALLEL_JOBS}" ]] || positive_int PARALLEL_JOBS "${PARALLEL_JOBS}"

# -----------------------------
# validate
# -----------------------------
[[ -d "$ASSEMBLY_DIR" ]] || die "ASSEMBLY_DIR not found: $ASSEMBLY_DIR"
[[ -d "$READS_DIR" ]]    || die "READS_DIR not found: $READS_DIR"
mkdir -p "$WORKDIR"

WORKDIR="$(realpath "$WORKDIR")"
ASSEMBLY_DIR="$(realpath "$ASSEMBLY_DIR")"
READS_DIR="$(realpath "$READS_DIR")"

# nextPolish check
if [[ ! -x "$NEXTPOLISH" ]]; then
  # allow conda-installed nextPolish in PATH
  if command -v "$NEXTPOLISH" >/dev/null 2>&1; then
    NEXTPOLISH="$(command -v "$NEXTPOLISH")"
  else
    die "nextPolish not executable: $NEXTPOLISH (use --nextpolish PATH)"
  fi
fi

NEXTPOLISH=$(realpath "$NEXTPOLISH")

# numeric checks
[[ "$THREADS" =~ ^[0-9]+$ ]] || die "--threads must be integer"
[[ "$PARALLEL_JOBS" =~ ^[0-9]+$ ]] || die "--jobs must be integer"
(( THREADS >= 1 )) || die "--threads must be >= 1"
(( PARALLEL_JOBS >= 1 )) || die "--jobs must be >= 1"

echo "=============================="
echo "NextPolish batch configuration"
echo "  assembly_dir   : $ASSEMBLY_DIR"
echo "  reads_dir      : $READS_DIR"
echo "  outdir         : $WORKDIR"
echo "  threads/job    : $THREADS"
echo "  parallel_jobs  : $PARALLEL_JOBS"
echo "  asm_glob       : $ASM_GLOB"
echo "  asm_suffix     : $ASM_SUFFIX"
echo "  r1_regex       : $R1_REGEX"
echo "  r2_regex       : $R2_REGEX"
echo "  nextpolish     : $NEXTPOLISH"
echo "  reads_mode     : $([[ "$COPY_READS" == true ]] && echo copy || echo symlink)"
echo "  keep_sampledir : $KEEP_SAMPLE_DIR"
echo "=============================="

# -----------------------------
# main loop
PIDS=()
# -----------------------------
shopt -s nullglob
assemblies=( "$ASSEMBLY_DIR"/$ASM_GLOB )
(( ${#assemblies[@]} > 0 )) || die "No assemblies matched: $ASSEMBLY_DIR/$ASM_GLOB"

for ASSEMBLY in "${assemblies[@]}"; do
  [[ -f "$ASSEMBLY" ]] || continue

  base="$(basename "$ASSEMBLY")"
  if [[ "$base" != *"$ASM_SUFFIX" ]]; then
    # if suffix doesn't match, fallback to filename without extension
    PREFIX="${base%.*}"
  else
    PREFIX="${base%$ASM_SUFFIX}"
  fi

  OUT_FASTA="$WORKDIR/${PREFIX}.nextpolish.fasta"
  if [[ -s "$OUT_FASTA" ]]; then
    echo "⏭️  $PREFIX already exists: $OUT_FASTA (skip)"
    continue
  fi

  # find reads
  R1=$(find_unique_read "$READS_DIR" "$PREFIX" "$R1_REGEX")
  R2=$(find_unique_read "$READS_DIR" "$PREFIX" "$R2_REGEX")

  if [[ -z "${R1:-}" || -z "${R2:-}" || ! -f "$R1" || ! -f "$R2" ]]; then
    echo "⚠️  Reads not found for $PREFIX, skipping..."
    echo "    looked for: .*/${PREFIX}.*${R1_REGEX}  and  .*/${PREFIX}.*${R2_REGEX}"
    continue
  fi

  SAMPLE_DIR="$WORKDIR/$PREFIX"
  mkdir -p "$SAMPLE_DIR"

  # put reads into sample dir (copy or symlink)
  if [[ "$COPY_READS" == true ]]; then
    cp -f "$R1" "$R2" "$SAMPLE_DIR/"
  else
    ln -sf "$R1" "$SAMPLE_DIR/"
    ln -sf "$R2" "$SAMPLE_DIR/"
  fi

  R1_BN="$(basename "$R1")"
  R2_BN="$(basename "$R2")"
  ASSEMBLY_ABS="$(realpath "$ASSEMBLY")"

  echo "🚀 Starting NextPolish for $PREFIX"
  echo "   Assembly: $ASSEMBLY_ABS"
  echo "   Reads:    $R1_BN, $R2_BN"
  echo "   Workdir:  $SAMPLE_DIR"

  # sgs.fofn
  FOFN="$SAMPLE_DIR/sgs.fofn"
  printf "%s\n%s\n" "$R1_BN" "$R2_BN" > "$FOFN"

  # run.cfg
  CFG="$SAMPLE_DIR/run.cfg"
  cat > "$CFG" <<EOF
[General]
job_type = local
job_prefix = nextPolish
task = best
rewrite = yes
deltmp = yes
rerun = 3
parallel_jobs = 2
multithread_jobs = $THREADS
genome = $ASSEMBLY_ABS
genome_size = auto
workdir = $SAMPLE_DIR
polish_options = -p {multithread_jobs}

[sgs_option]
sgs_fofn = $FOFN
sgs_options = -max_depth 100 -bwa
EOF

  (
    set -euo pipefail
    cd "$SAMPLE_DIR"

    "$NEXTPOLISH" "$CFG"

    FINAL_FASTA="$(find "$SAMPLE_DIR" -type f -name "genome.nextpolish.fasta" | head -n 1 || true)"
    if [[ -f "$FINAL_FASTA" ]]; then
      cp -f "$FINAL_FASTA" "$WORKDIR/${PREFIX}.nextpolish.fasta"
      echo "✅ Finished: $WORKDIR/${PREFIX}.nextpolish.fasta"

      if [[ "$KEEP_SAMPLE_DIR" == false ]]; then
        rm -rf "$SAMPLE_DIR"
      fi
    else
      echo "❌ No polished fasta found for $PREFIX"
      exit 1
    fi
  ) &
  PIDS+=("$!")

  # concurrency control
  while (( $(jobs -r | wc -l) >= PARALLEL_JOBS )); do
    sleep 5
  done
done

FAILED=0
for pid in "${PIDS[@]}"; do
  wait "$pid" || FAILED=1
done
(( FAILED == 0 )) || { echo "ERROR: One or more samples failed" >&2; exit 1; }
echo "🎉 All polishing jobs finished!"
