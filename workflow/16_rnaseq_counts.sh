#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# STAR + featureCounts pipeline (multi-species)
# Input: a TXT/TSV with 5 columns (TAB-separated recommended)
#   1) species
#   2) R1 fastq(.gz)
#   3) R2 fastq(.gz)
#   4) genome fasta
#   5) annotation gtf
#
# Output:
#   OUTDIR/counts/{species}_count.txt   (gene x sample counts)
# plus:
#   OUTDIR/align/{species}/*.sorted.bam
#   OUTDIR/star_index/{species}/
#   OUTDIR/logs/
#
# Requirements in PATH: STAR, samtools, featureCounts
# ============================================================

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/lib/common.sh"
if [[ "${1:-}" == -h || "${1:-}" == --help ]]; then
  cat <<'EOF'
Usage: 16_rnaseq_counts.sh MANIFEST.tsv OUTDIR
Manifest: five TAB-separated columns, with a # comment header:
  species  R1  R2  genome_fasta  annotation_gtf
Paths resolve relative to the launch directory. Each species needs one genome/GTF pair.
Environment: THREADS=24 FC_THREADS=12 SJDB_OVERHANG=149 STRAND=0
            FC_FEATURE_TYPE=exon FC_GENE_ATTR=gene_id COUNT_MODE=legacy
COUNT_MODE: legacy (original command), reads (-p), fragments (-p --countReadPairs).
COUNT_MODE can be legacy, reads, or fragments.
Existing indices/BAMs are reused; use a fresh output directory after input changes.
EOF
  exit 0
fi
(( $# == 2 )) || { echo "ERROR: Expected MANIFEST.tsv OUTDIR (see --help)" >&2; exit 2; }
# --------- user settings ----------
INPUT_TXT="${1:-meta/runs.txt}"     # your 5-column file
OUTDIR="${2:-results_star_fc}"      # output dir

THREADS="${THREADS:-24}"            # STAR threads per sample
FC_THREADS="${FC_THREADS:-12}"      # featureCounts threads
SJDB_OVERHANG="${SJDB_OVERHANG:-149}"  # readLength-1 (PE150 -> 149, PE100 -> 99)
STRAND="${STRAND:-0}"               # 0=unstranded, 1=stranded, 2=reverse stranded
FC_FEATURE_TYPE="${FC_FEATURE_TYPE:-exon}"
FC_GENE_ATTR="${FC_GENE_ATTR:-gene_id}"

# STAR alignment tuning for yeast-like genomes
STAR_EXTRA_ALIGN_OPTS=(
  --outFilterMultimapNmax 10
  --outFilterMismatchNoverLmax 0.04
  --alignIntronMin 20
  --alignIntronMax 2000
  --alignMatesGapMax 2000
  --outSAMattributes NH HI AS nM MD
)

COUNT_MODE="${COUNT_MODE:-legacy}"
FC_PAIR_OPTS=()
case "$COUNT_MODE" in
  legacy) ;;
  reads) FC_PAIR_OPTS=(-p);;
  fragments) FC_PAIR_OPTS=(-p --countReadPairs);;
  *) echo "ERROR: COUNT_MODE must be legacy, reads, or fragments" >&2; exit 2;;
esac
positive_int THREADS "$THREADS"
positive_int FC_THREADS "$FC_THREADS"
[[ "$STRAND" =~ ^[012]$ ]] || { echo "ERROR: STRAND must be 0, 1, or 2" >&2; exit 2; }
python3 "$SCRIPT_DIR/lib/rna_io.py" validate "$INPUT_TXT"

# --------- folders ----------
IDXDIR="${OUTDIR}/star_index"
ALNDIR="${OUTDIR}/align"
COUNTDIR="${OUTDIR}/counts"
LOGDIR="${OUTDIR}/logs"

mkdir -p "$OUTDIR" "$IDXDIR" "$ALNDIR" "$COUNTDIR" "$LOGDIR"

# --------- checks ----------
need_cmd() { command -v "$1" >/dev/null 2>&1 || { echo "[ERROR] missing command: $1" >&2; exit 1; }; }
need_file() { [[ -s "$1" ]] || { echo "[ERROR] missing/empty file: $1" >&2; exit 1; }; }

need_cmd STAR
need_cmd samtools
need_cmd featureCounts
need_file "$INPUT_TXT"

# Determine if a file is gzipped by suffix
is_gz() { [[ "$1" == *.gz ]]; }

# Generate a sample name from R1 path (strip directory + common R1 suffixes)
sample_from_r1() {
  local r1="$1"
  local b
  b="$(basename "$r1")"
  b="${b%.fastq.gz}"; b="${b%.fq.gz}"; b="${b%.fastq}"; b="${b%.fq}"
  # strip common read1 markers (best-effort)
  b="${b%_R1}"; b="${b%_1}"; b="${b%R1}"; b="${b%1}"
  echo "$b"
}

# --------- parse unique species ----------
mapfile -t SPECIES_LIST < <(awk '
  BEGIN{FS="\t"}
  $0 ~ /^[[:space:]]*$/ {next}
  $0 ~ /^[[:space:]]*#/ {next}
  {print $1}
' "$INPUT_TXT" | sort -u)

if [[ ${#SPECIES_LIST[@]} -eq 0 ]]; then
  echo "[ERROR] No species found in $INPUT_TXT. Ensure it is TAB-separated with 5 columns." >&2
  exit 1
fi

echo "[INFO] Found species (${#SPECIES_LIST[@]}): ${SPECIES_LIST[*]}"
echo "[INFO] Input: $INPUT_TXT"
echo "[INFO] Output: $OUTDIR"

# --------- build STAR index per species ----------
build_index_for_species() {
  local sp="$1"

  # pick the first line for this species to get genome/gtf paths
  local fa gtf
  fa="$(awk -v sp="$sp" 'BEGIN{FS="\t"}
      $0 ~ /^[[:space:]]*$/ {next}
      $0 ~ /^[[:space:]]*#/ {next}
      $1==sp {print $4; exit}
    ' "$INPUT_TXT")"
  gtf="$(awk -v sp="$sp" 'BEGIN{FS="\t"}
      $0 ~ /^[[:space:]]*$/ {next}
      $0 ~ /^[[:space:]]*#/ {next}
      $1==sp {print $5; exit}
    ' "$INPUT_TXT")"

  [[ -n "$fa" && -n "$gtf" ]] || { echo "[ERROR] Cannot find genome/gtf for species=$sp in $INPUT_TXT" >&2; exit 1; }
  need_file "$fa"
  need_file "$gtf"

  local idx="${IDXDIR}/${sp}"
  mkdir -p "$idx"

  if [[ -s "${idx}/Genome" ]]; then
    echo "[INFO] STAR index exists: $sp -> $idx"
    return 0
  fi

  echo "[INFO] Building STAR index: $sp"
  STAR \
    --runThreadN "$THREADS" \
    --runMode genomeGenerate \
    --genomeDir "$idx" \
    --genomeFastaFiles "$fa" \
    --sjdbGTFfile "$gtf" \
    --sjdbOverhang "$SJDB_OVERHANG" \
    > "${LOGDIR}/STAR_genomeGenerate.${sp}.log" 2>&1
}

# --------- align one row (one sample) ----------
align_one_row() {
  local sp="$1" r1="$2" r2="$3"
  need_file "$r1"
  need_file "$r2"

  local idx="${IDXDIR}/${sp}"
  [[ -s "${idx}/Genome" ]] || { echo "[ERROR] STAR index missing for $sp at $idx" >&2; exit 1; }

  local sample
  sample="$(sample_from_r1 "$r1")"

  local out_sp="${ALNDIR}/${sp}"
  mkdir -p "$out_sp"

  local bam="${out_sp}/${sample}.sorted.bam"
  local prefix="${out_sp}/${sample}."

  if [[ -s "$bam" && -s "${bam}.bai" ]]; then
    echo "[INFO] Skip align (exists): $sp / $sample"
    return 0
  fi

  echo "[INFO] Align: $sp / $sample"
  local readcmd=()
  if is_gz "$r1"; then readcmd=(--readFilesCommand zcat); fi

  STAR \
    --runThreadN "$THREADS" \
    --genomeDir "$idx" \
    --readFilesIn "$r1" "$r2" \
    "${readcmd[@]}" \
    --outFileNamePrefix "$prefix" \
    --outSAMtype BAM SortedByCoordinate \
    --limitBAMsortRAM 16000000000 \
    "${STAR_EXTRA_ALIGN_OPTS[@]}" \
    > "${LOGDIR}/STAR_align.${sp}.${sample}.log" 2>&1

  mv -f "${prefix}Aligned.sortedByCoord.out.bam" "$bam"
  samtools index -@ "$THREADS" "$bam"

  # keep mapping summary
  [[ -s "${prefix}Log.final.out" ]] && cp -f "${prefix}Log.final.out" "${LOGDIR}/STAR_final.${sp}.${sample}.txt" || true

  # cleanup noisy intermediates
  rm -f "${prefix}Aligned.out.sam" "${prefix}SJ.out.tab" "${prefix}Log.out" "${prefix}Log.progress.out" 2>/dev/null || true
}

# --------- run featureCounts per species ----------
count_species() {
  local sp="$1"

  # genome/gtf from first line of this species
  local gtf
  gtf="$(awk -v sp="$sp" 'BEGIN{FS="\t"}
      $0 ~ /^[[:space:]]*$/ {next}
      $0 ~ /^[[:space:]]*#/ {next}
      $1==sp {print $5; exit}
    ' "$INPUT_TXT")"
  [[ -n "$gtf" ]] || { echo "[ERROR] Cannot find gtf for $sp in $INPUT_TXT" >&2; exit 1; }
  need_file "$gtf"

  local bam_dir="${ALNDIR}/${sp}"
  [[ -d "$bam_dir" ]] || { echo "[WARN] No bam dir for $sp; skip counts"; return 0; }

  # Count only samples in this manifest; exclude stale BAMs from previous runs.
  local bams=() manifest_r1 sample
  while IFS= read -r manifest_r1; do
    sample=$(sample_from_r1 "$manifest_r1")
    bams+=("${bam_dir}/${sample}.sorted.bam")
    need_file "${bams[-1]}"
  done < <(awk -F '\t' -v sp="$sp" '$1==sp {print $2}' "$INPUT_TXT")
  if [[ ${#bams[@]} -eq 0 ]]; then
    echo "[WARN] No BAMs for $sp; skip counts"
    return 0
  fi

  local fc_raw="${COUNTDIR}/${sp}.featureCounts.txt"
  local fc_sum="${COUNTDIR}/${sp}.featureCounts.summary.txt"
  local out_mat="${COUNTDIR}/${sp}_count.txt"

  echo "[INFO] featureCounts: $sp (n_bams=${#bams[@]})"
  featureCounts \
    "${FC_PAIR_OPTS[@]}" \
    -T "$FC_THREADS" \
    -s "$STRAND" \
    -a "$gtf" \
    -o "$fc_raw" \
    -t "$FC_FEATURE_TYPE" \
    -g "$FC_GENE_ATTR" \
    "${bams[@]}" \
    > "${LOGDIR}/featureCounts.${sp}.log" 2>&1

  [[ -s "${fc_raw}.summary" ]] && mv -f "${fc_raw}.summary" "$fc_sum" || true

  # Clean matrix: Geneid + sample columns, strip bam path & suffix.
  python3 "$SCRIPT_DIR/lib/rna_io.py" counts "$fc_raw" "$out_mat"

  echo "[INFO] Wrote counts: $out_mat"
}

# ============================================================
# MAIN
# ============================================================

echo "[INFO] Step 1/3: Build indices"
for sp in "${SPECIES_LIST[@]}"; do
  build_index_for_species "$sp"
done

echo "[INFO] Step 2/3: Align samples"
# Read input rows: species r1 r2 fasta gtf
# (We ignore fasta/gtf per row after index is built, but they must be consistent.)
awk '
  BEGIN{FS="\t"}
  $0 ~ /^[[:space:]]*$/ {next}
  $0 ~ /^[[:space:]]*#/ {next}
  {print $1 "\t" $2 "\t" $3}
' "$INPUT_TXT" | while IFS=$'\t' read -r sp r1 r2; do
  align_one_row "$sp" "$r1" "$r2"
done

echo "[INFO] Step 3/3: featureCounts per species"
for sp in "${SPECIES_LIST[@]}"; do
  count_species "$sp"
done

echo "[INFO] DONE."
echo "  Counts in: ${COUNTDIR}  (files: {species}_count.txt)"
