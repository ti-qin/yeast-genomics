#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/lib/common.sh"
DATA_DIR=""
OUTDIR=""
COMBINE_RM="${COMBINE_RM:-/opt/RepeatMasker/util/combineRMFiles.pl}"
RM1_THREADS=128
RM2_THREADS=24
RM_THREADS=32
usage() {
  cat <<'EOF'
Usage: 12_repeat_annotation.sh -i GENOME_DIR -o OUTDIR [options]
Input: top-level *.nextpolish.fasta files.
  --modeler-threads INT    RepeatModeler threads (default: 128)
  --denovo-threads INT     De novo RepeatMasker -pa value (default: 24)
  --species-threads INT    Species RepeatMasker -pa value (default: 32)
  --combine-rm FILE        combineRMFiles.pl (default: /opt/RepeatMasker/util/combineRMFiles.pl)
  -h, --help              Show help
Run inside a configured TETools environment with repeat libraries installed.
EOF
}
validate_value_options "-i -o --modeler-threads --denovo-threads --species-threads --combine-rm" "$@"
while (( $# )); do
  case "$1" in
    -i) DATA_DIR="$2"; shift 2;;
    -o) OUTDIR="$2"; shift 2;;
    --modeler-threads) RM1_THREADS="$2"; shift 2;;
    --denovo-threads) RM2_THREADS="$2"; shift 2;;
    --species-threads) RM_THREADS="$2"; shift 2;;
    --combine-rm) COMBINE_RM="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) usage >&2; exit 2;;
  esac
done
[[ -d "$DATA_DIR" && -n "$OUTDIR" ]] || { usage >&2; exit 2; }
positive_int RM1_THREADS "$RM1_THREADS"
positive_int RM2_THREADS "$RM2_THREADS"
positive_int RM_THREADS "$RM_THREADS"
for cmd in BuildDatabase RepeatModeler RepeatMasker perl rmOutToGFF3.pl maskFile.pl; do
  command -v "$cmd" >/dev/null || { echo "ERROR: Missing $cmd" >&2; exit 1; }
done
[[ -f "$COMBINE_RM" ]] || { echo "ERROR: Missing $COMBINE_RM" >&2; exit 1; }
COMBINE_RM=$(realpath "$COMBINE_RM")
mkdir -p "$OUTDIR"
DATA_DIR=$(realpath "$DATA_DIR")
OUTDIR=$(realpath "$OUTDIR")
shopt -s nullglob
FASTAS=( "$DATA_DIR"/*.nextpolish.fasta )
(( ${#FASTAS[@]} > 0 )) || { echo "ERROR: No *.nextpolish.fasta inputs" >&2; exit 1; }

# 遍历所有 *.nextpolish.fasta
for fasta in "${FASTAS[@]}"; do
    [[ -e "$fasta" ]] || continue  # no files

    # 取出样本名 G8699
    sample=$(basename "$fasta" .nextpolish.fasta)
    echo "==== Processing sample: $sample ===="

    # 建立样本工作目录
    SAMPLE_DIR="$OUTDIR/$sample"
    mkdir -p "$SAMPLE_DIR"
    cd "$SAMPLE_DIR"

    # 软链接 fasta（避免复制大文件）
    ln -sf "$fasta" "./${sample}.nextpolish.fasta"

    echo "[1] BuildDatabase"
    BuildDatabase -name "$sample" "./${sample}.nextpolish.fasta"

    echo "[2] RepeatModeler"
    RepeatModeler -database "$sample" -threads $RM1_THREADS -LTRStruct \
        1>repeatmodeler.log 2>repeatmodeler.err

    # RepeatMasker Step 1 (species)
    mkdir -p repeatmasker_dir
    echo "[3] RepeatMasker (species saccharomycotina)"
    RepeatMasker -e rmblast -pa $RM_THREADS -gff -q -a \
        -species saccharomycotina \
        -dir ./repeatmasker_dir \
        "./${sample}.nextpolish.fasta" \
        1>repeatmasker_spe.log 2>repeatmasker_spe.err

    # RepeatMasker Step 2 (denovo)
    mkdir -p repeatmodel_dir
    echo "[4] RepeatMasker (denovo lib)"
    RepeatMasker -e rmblast -pa $RM2_THREADS -gff -q -a \
        -lib "${sample}-families.fa" \
        -dir ./repeatmodel_dir \
        "./repeatmasker_dir/${sample}.nextpolish.fasta.masked" \
        1>repeatmasker_denovo.log 2>repeatmasker_denovo.err

    # Combine
    echo "[5] combineRMFiles"
    mkdir -p repeatmaskerout_combine
    perl "$COMBINE_RM" \
        "./repeatmasker_dir/${sample}.nextpolish.fasta" \
        "./repeatmodel_dir/${sample}.nextpolish.fasta.masked" \
        "repeatmaskerout_combine/${sample}"

    # Convert OUT → GFF3
    echo "[6] rmOutToGFF3"
    rmOutToGFF3.pl "repeatmaskerout_combine/${sample}.out" \
        > "repeatmaskerout_combine/${sample}.out.gff"

    # Filter
    echo "[7] Filter OUT"
    awk '!/Low_complexity|Simple_repeat|Unknown/' \
      "repeatmaskerout_combine/${sample}.out" \
      > "repeatmaskerout_combine/${sample}.filter.out"

    # Final softmasking
    echo "[8] Softmask final FASTA"
    maskFile.pl -fasta "./${sample}.nextpolish.fasta" \
        -annotations "repeatmaskerout_combine/${sample}.filter.out" \
        -softmask

    echo "==== $sample DONE ===="
    cd - >/dev/null
done

echo "All samples finished."