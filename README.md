# Yeast genome assembly and annotation workflow

Code used for genome assembly, polishing, quality assessment, annotation, and RNA-seq quantification in:

> **400 million years of evolution reveal the determinants of gene order conservation**

Thirty budding yeast genomes were assembled independently with Canu, NECAT, and NextDenovo. Draft assemblies were filtered, polished with long and short reads, assessed, and one assembly was selected for each species before annotation.

## Repository structure

```text
workflow/
├── 01_assembly_canu.sh
├── 02_assembly_necat.sh
├── 03_assembly_nextdenovo.sh
├── 04_remove_redundant_contigs.sh
├── 05_remove_mitochondrial_contigs.sh
├── 06_assembly_statistics.sh
├── 07_polish_racon.sh
├── 08_polish_nextpolish.sh
├── 09_assess_merqury.sh
├── 10_assess_compleasm.sh
├── 11_reference_alignment.sh
├── 12_repeat_annotation.sh
├── 13_gene_annotation.sh
├── 14_assess_annotation.sh
├── 15_mash_distance.sh
├── 16_rnaseq_counts.sh
├── examples/
├── lib/
└── references/
```

The scripts are numbered in execution order. Assembly statistics and assessment scripts are used to compare the three assembly branches before selecting the final assembly.

## Workflow

```text
Long reads
   ├── Canu
   ├── NECAT
   └── NextDenovo
          ↓
Redundant-contig removal
          ↓
Mitochondrial-contig removal
          ↓
Racon polishing
          ↓
NextPolish polishing
          ↓
Assembly assessment and selection
          ↓
RepeatMasker / RepeatModeler
          ↓
BRAKER3 annotation
          ↓
Annotation assessment and RNA-seq quantification
```

Redundant contigs are removed at 80% coverage by other contigs. Putative mitochondrial contigs are removed when at least two of the following criteria are met: high sequencing coverage (MAD-based score > 5), a mitochondrial nucleotide alignment of at least 500 bp and 80% identity, or matches to at least three mitochondrial proteins.

## Usage

Run scripts from the repository root. Each script provides its complete command-line interface:

```bash
workflow/01_assembly_canu.sh --help
workflow/02_assembly_necat.sh --help
workflow/03_assembly_nextdenovo.sh --help
```

Example:

```bash
workflow/03_assembly_nextdenovo.sh \
  -r rawdata/long_reads \
  -w results/01_assembly/nextdenovo \
  -g 12000000 \
  -t 16 \
  -j 2
```

Input sequencing reads, BUSCO databases, repeat libraries, and container images are not included. The S288C mitochondrial reference and six protein markers used in step 05 are provided in `workflow/references/`. File paths and computing resources can be adjusted according to the actual server environment. RNA-seq input is specified with the template in `workflow/examples/rna_manifest.tsv`.

## Main dependencies

Canu, NECAT, NextDenovo, BLAST+, minimap2, samtools, SeqKit, MUMmer4, Racon, NextPolish, Merqury, Compleasm, RepeatModeler, RepeatMasker, BRAKER3, Mash, STAR, featureCounts, Python 3, pandas, NumPy, and Biopython.

## Attribution

The redundancy-removal utilities in `workflow/lib/` were adapted from [HaploTeam/1086YeastGenomes](https://github.com/HaploTeam/1086YeastGenomes), associated with *From genotype to phenotype with 1,086 near telomere-to-telomere yeast genomes*.
