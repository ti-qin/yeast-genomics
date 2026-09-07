#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/lib/common.sh"
validate_value_options "--fasta-dir --fastq-dir --keep-tmp --len-max --mad-z --mt-aa --mt-dna --nuc-min-idy --nuc-min-len --outdir --tbl-evalue --tbl-min-bits --threads -a -d -i -o -q -t" "$@"

# -------------------------
# defaults
# -------------------------
FASTA_DIR=""
FASTQ_DIR=""
OUT_DIR="./results/clean/nuclear"
MT_DNA="$SCRIPT_DIR/references/S288C_mitochondrion_NC_001224.1.fasta"
MT_AA="$SCRIPT_DIR/references/S288C_mitochondrial_markers_6.faa"
THREADS=32

LEN_MAX=200000

# nucmer counting thresholds
NUC_MIN_LEN=500
NUC_MIN_IDY=80

# tblastn thresholds
TBL_EVALUE="1e-5"
TBL_MIN_BITSCORE=50

# MAD threshold (modified z-score)
MAD_Z=5

# keep temporary directories
KEEP_TMP=true

usage() {
cat <<EOF
Usage: $(basename "$0") -i FASTA_DIR -q FASTQ_DIR [options]

Required:
  -i, --fasta-dir DIR     Directory containing *NR.fasta (e.g., G162.flye.NR.fasta)
  -q, --fastq-dir DIR     Directory containing reads (e.g., G162.fastq.gz)
Optional:
  -d, --mt-dna FILE       Mitochondrial DNA reference (default: bundled S288C sequence)
  -a, --mt-aa FILE        Mitochondrial proteins (default: bundled S288C markers)
  -o, --outdir DIR        Output folder (default: ./results/clean/nuclear)
  -t, --threads INT       Threads (default: 32)
      --len-max INT       Only run nucmer/tblastn on contigs <= this (default: 200000)
      --nuc-min-len INT   Count nucmer blocks with aligned LEN2 >= this (default: 500)
      --nuc-min-idy INT   Count nucmer blocks with IDY >= this (default: 80)
      --tbl-evalue STR    tblastn evalue (default: 1e-5)
      --tbl-min-bits INT  tblastn min bitscore (default: 50)
      --mad-z FLOAT       MAD modified z-score threshold for high-cov outliers (default: 5)
      --keep-tmp BOOL     Keep temp dirs true|false (default: true)

Classification:
  For contigs <= --len-max (default: 200000 bp), require at least TWO:
  high coverage modified-z > --mad-z (default: 5), at least one accepted
  nucleotide block, or at least three distinct mitochondrial protein hits.

Outputs (always in OUTDIR root):
  <outdir>/mito_screen.summary.tsv
  <outdir>/<FASTA_BASENAME_NOEXT>.nuclear.fasta
  <outdir>/<FASTA_BASENAME_NOEXT>.mito_contigs.list
EOF
}

# -------------------------
# parse args
# -------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--fasta-dir) FASTA_DIR="$2"; shift 2;;
    -q|--fastq-dir) FASTQ_DIR="$2"; shift 2;;
    -d|--mt-dna) MT_DNA="$2"; shift 2;;
    -a|--mt-aa) MT_AA="$2"; shift 2;;
    -o|--outdir) OUT_DIR="$2"; shift 2;;
    -t|--threads) THREADS="$2"; shift 2;;
    --len-max) LEN_MAX="$2"; shift 2;;
    --nuc-min-len) NUC_MIN_LEN="$2"; shift 2;;
    --nuc-min-idy) NUC_MIN_IDY="$2"; shift 2;;
    --tbl-evalue) TBL_EVALUE="$2"; shift 2;;
    --tbl-min-bits) TBL_MIN_BITSCORE="$2"; shift 2;;
    --mad-z) MAD_Z="$2"; shift 2;;
    --keep-tmp) KEEP_TMP="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "ERROR: Unknown option: $1"; usage; exit 1;;
  esac
done
[[ -z "${THREADS}" ]] || positive_int THREADS "${THREADS}"

# -------------------------
# checks
# -------------------------
[[ -n "$FASTA_DIR" ]] || { echo "ERROR: --fasta-dir required"; usage; exit 1; }
[[ -n "$FASTQ_DIR" ]] || { echo "ERROR: --fastq-dir required"; usage; exit 1; }
[[ -f "$MT_DNA" ]] || { echo "ERROR: mt_dna not found: $MT_DNA"; exit 1; }
[[ -f "$MT_AA"  ]] || { echo "ERROR: mt_aa not found: $MT_AA"; exit 1; }

for cmd in minimap2 samtools seqkit nucmer delta-filter show-coords makeblastdb tblastn python3; do
  command -v "$cmd" >/dev/null 2>&1 || { echo "ERROR: $cmd not in PATH"; exit 1; }
done

if [[ "$KEEP_TMP" != "true" && "$KEEP_TMP" != "false" ]]; then
  echo "ERROR: --keep-tmp must be true or false"
  exit 1
fi

mkdir -p "$OUT_DIR"

SUMMARY="$OUT_DIR/mito_screen.summary.tsv"
echo -e "fasta\tcontig\tcov\tmean_cov\tlen\tnucmer_blocks\ttblastn_genes\tcov_outlier\tmito_flag" > "$SUMMARY"

# -------------------------
# helpers
# -------------------------
find_reads() {
  local sample="$1"
  local fq=""
  for ext in fastq.gz fq.gz fastq fq; do
    if [[ -f "$FASTQ_DIR/${sample}.${ext}" ]]; then
      fq="$FASTQ_DIR/${sample}.${ext}"
      break
    fi
  done
  echo "$fq"
}

mad_flag_cov_outlier() {
  # input: contig\tlen\tcov
  # output: contig\tlen\tcov\tmean_cov\tcov_outlier
  local in_tsv="$1"
  local out_tsv="$2"
  local mad_z="$3"

  python3 - "$in_tsv" "$out_tsv" "$mad_z" <<'PY'
import sys
import numpy as np
import pandas as pd

in_tsv, out_tsv, mad_z = sys.argv[1], sys.argv[2], float(sys.argv[3])

df = pd.read_csv(in_tsv, sep="\t", header=None, names=["contig","len","cov"])
if df.empty:
    raise SystemExit("EMPTY_CONTIG_TABLE")

df["cov"] = df["cov"].astype(float)
mean_cov = float(df["cov"].mean())

cov = df["cov"].to_numpy()
med = float(np.median(cov))
mad = float(np.median(np.abs(cov - med)))

if mad == 0.0:
    outlier = np.zeros_like(cov, dtype=int)
else:
    mod_z = 0.6745 * (cov - med) / mad
    outlier = (mod_z > mad_z).astype(int)  # only high outliers

df["mean_cov"] = mean_cov
df["cov_outlier"] = outlier

df.to_csv(out_tsv, sep="\t", index=False, header=False)
PY
}

# -------------------------
# main
# -------------------------
shopt -s nullglob
FA_FILES=("$FASTA_DIR"/*NR.fasta)
shopt -u nullglob
[[ ${#FA_FILES[@]} -gt 0 ]] || { echo "ERROR: No *NR.fasta in $FASTA_DIR"; exit 1; }

echo "=== Start mito screen (map-based cov + MAD): $(date) ==="
echo "FASTA_DIR: $FASTA_DIR"
echo "FASTQ_DIR: $FASTQ_DIR"
echo "OUT_DIR:   $OUT_DIR"
echo "threads:   $THREADS"
echo "len_max:   $LEN_MAX"
echo "MAD_Z:     $MAD_Z"
echo "keep_tmp:  $KEEP_TMP"
echo "==============================================="

for fa in "${FA_FILES[@]}"; do
  base=$(basename "$fa")              # G162.flye.NR.fasta
  name="${base%.fasta}"               # G162.flye.NR
  sample="${name%%.*}"                # G162

  reads=$(find_reads "$sample")
  if [[ -z "$reads" ]]; then
    echo "❌ ERROR: reads not found for sample=$sample (expected $FASTQ_DIR/${sample}.fastq.gz etc.)"
    exit 1
  fi

  echo ">>> $base"
  echo "    sample: $sample"
  echo "    reads : $reads"

  # temp dir per fasta
  sdir="$OUT_DIR/$name"
  mkdir -p "$sdir"

  # final outputs in OUT_DIR root
  nuclear_out="$OUT_DIR/${name}.nuclear.fasta"
  mito_list_out="$OUT_DIR/${name}.mito_contigs.list"

  # contig lengths
  lenmap="$sdir/contig.len.tsv"
  seqkit fx2tab -n -l "$fa" | awk 'BEGIN{OFS="\t"}{print $1,$NF}' > "$lenmap"

  # map reads -> bam
  bam="$sdir/reads_to_asm.bam"
  minimap2 -t "$THREADS" -ax map-ont "$fa" "$reads" \
    | samtools sort -@ "$THREADS" -o "$bam" -
  samtools index "$bam"

  # mean coverage per contig
  cov_tsv="$sdir/contig.cov.tsv"  # contig\tcov
  samtools depth -aa "$bam" \
    | awk 'BEGIN{FS=OFS="\t"} {sum[$1]+=$3; n[$1]++} END{for (c in n) print c, sum[c]/n[c]}' \
    | sort -k1,1 > "$cov_tsv"

  # merge len + cov => contig len cov
  covlen_tsv="$sdir/contig.len.cov.tsv"
  join -t $'\t' -a 1 -e "0" -o 0,1.2,2.2 \
    <(sort -k1,1 "$lenmap") \
    <(sort -k1,1 "$cov_tsv") > "$covlen_tsv"

  # MAD outlier => contig len cov mean_cov outlier
  stats_tsv="$sdir/contig.stats.tsv"
  mad_flag_cov_outlier "$covlen_tsv" "$stats_tsv" "$MAD_Z"
  sort -k1,1 "$stats_tsv" -o "$stats_tsv"

  # short contigs list (<= LEN_MAX)
  short_list="$sdir/len_le_${LEN_MAX}.list"
  short_fa="$sdir/len_le_${LEN_MAX}.fa"
  awk -v m="$LEN_MAX" '$2<=m{print $1}' "$stats_tsv" > "$short_list"

  # If no short contigs: keep original fasta, and summary has nuc/tbl=0, mito_flag=0
  if [[ ! -s "$short_list" ]]; then
    echo "    (No contigs <= ${LEN_MAX} bp; keep original fasta as nuclear)"
    cp -f "$fa" "$nuclear_out"
    : > "$mito_list_out"

    # stats_tsv cols: contig len cov mean outlier
    while IFS=$'\t' read -r contig len cov mean outlier; do
      echo -e "${base}\t${contig}\t${cov}\t${mean}\t${len}\t0\t0\t${outlier}\t0" >> "$SUMMARY"
    done < "$stats_tsv"

    echo "    output: $nuclear_out"
    [[ "$KEEP_TMP" == "false" ]] && rm -rf "$sdir"
    continue
  fi

  seqkit grep -f "$short_list" "$fa" > "$short_fa"

  # 2) nucmer vs mt_dna => nucmer_blocks per contig
  nuc_counts="$sdir/nucmer.blocks.count.tsv"
  : > "$nuc_counts"

  nuc_prefix="$sdir/nucmer_vs_mt"
  nuc_delta="${nuc_prefix}.delta"
  nuc_coords="$sdir/nucmer.coords.tsv"

  nucmer --maxmatch --nosimplify -t "$THREADS" -p "$nuc_prefix" "$MT_DNA" "$short_fa" > /dev/null 2>&1
  if [[ -s "$nuc_delta" ]]; then
    delta-filter -1 "$nuc_delta" > "$sdir/nucmer.filtered.delta"
    show-coords -THrd "$sdir/nucmer.filtered.delta" > "$nuc_coords"

    awk -v minlen="$NUC_MIN_LEN" -v minidy="$NUC_MIN_IDY" '
      BEGIN{OFS="\t"}
      NF>=11{
        len2=$6; idy=$7; qry=$11;
        if (len2+0 >= minlen && idy+0 >= minidy) cnt[qry]++
      }
      END{ for (q in cnt) print q, cnt[q] }
    ' "$nuc_coords" | sort -k1,1 > "$nuc_counts"
  fi

  # 3) tblastn vs mt_aa => unique genes per contig
  tbl_counts="$sdir/tblastn.genes.count.tsv"
  : > "$tbl_counts"

  makeblastdb -in "$short_fa" -dbtype nucl > /dev/null 2>&1

  tbl_out="$sdir/tblastn.tsv"
  tblastn -query "$MT_AA" -db "$short_fa" \
    -evalue "$TBL_EVALUE" -num_threads "$THREADS" \
    -outfmt '6 qseqid sseqid evalue bitscore pident length' \
    -max_target_seqs 100000 > "$tbl_out"

  awk -v minbs="$TBL_MIN_BITSCORE" '
    BEGIN{OFS="\t"}
    NF>=6{
      q=$1; s=$2; bs=$4;
      if (bs+0 >= minbs) hit[s SUBSEP q]=1
    }
    END{
      for (k in hit){
        split(k,a,SUBSEP); s=a[1]; cnt[s]++
      }
      for (s in cnt) print s, cnt[s]
    }
  ' "$tbl_out" | sort -k1,1 > "$tbl_counts"

  # merge: stats + nuc + tbl
  # stats_tsv: contig len cov mean outlier
  merged="$sdir/metrics.merged.tsv"

  join -t $'\t' -a 1 -e "0" -o 0,1.2,1.3,1.4,1.5,2.2 \
    <(sort -k1,1 "$stats_tsv") <(sort -k1,1 "$nuc_counts") > "$sdir/_stats_nuc.tsv"
  # _stats_nuc: contig len cov mean outlier nuc

  join -t $'\t' -a 1 -e "0" -o 0,1.2,1.3,1.4,1.5,1.6,2.2 \
    <(sort -k1,1 "$sdir/_stats_nuc.tsv") <(sort -k1,1 "$tbl_counts") > "$merged"
  # merged: contig len cov mean outlier nuc tbl

  # apply deletion rules (only len<=LEN_MAX)
  final="$sdir/final.withflag.tsv"
  awk -v m="$LEN_MAX" -f "$SCRIPT_DIR/lib/classify_mito.awk" "$merged" > "$final"
  # final: contig cov mean len nuc tbl outlier mito

  # lists
  keep_list="$sdir/keep_contigs.list"
  awk -F'\t' '$8==0{print $1}' "$final" | sort -u > "$keep_list"
  awk -F'\t' '$8==1{print $1}' "$final" | sort -u > "$mito_list_out"

  # output nuclear fasta
  seqkit grep -f "$keep_list" "$fa" > "$nuclear_out"

  # append to global summary (has cov & mean!)
  while IFS=$'\t' read -r contig cov mean len nuc tbl outlier mito; do
    echo -e "${base}\t${contig}\t${cov}\t${mean}\t${len}\t${nuc}\t${tbl}\t${outlier}\t${mito}" >> "$SUMMARY"
  done < "$final"

  echo "    mito contigs: $(wc -l < "$mito_list_out" | tr -d ' ')"
  echo "    output: $nuclear_out"

  [[ "$KEEP_TMP" == "false" ]] && rm -rf "$sdir"
done

echo "=== DONE: $(date) ==="
echo "Summary: $SUMMARY"
