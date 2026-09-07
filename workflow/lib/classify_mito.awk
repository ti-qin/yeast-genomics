# Input: contig length coverage mean_coverage coverage_outlier nuc_blocks protein_hits
# Output: contig coverage mean_coverage length nuc_blocks protein_hits outlier mito_flag
# Manuscript rule confirmed by the author: require at least two evidence classes.
BEGIN { FS=OFS="\t" }
{
    contig=$1; len=$2+0; cov=$3; mean=$4; out=$5+0; nuc=$6+0; tbl=$7+0;
    evidence=(out==1) + (nuc>=1) + (tbl>=3);
    mito=(len<=m && evidence>=2) ? 1 : 0;
    print contig, cov, mean, len, nuc, tbl, out, mito;
}
