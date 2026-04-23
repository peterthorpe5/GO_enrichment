# Barley GO enrichment for custom BaRT annotations

This repository contains an R script for running Gene Ontology enrichment on custom barley annotations, including BaRT-style transcript and gene identifiers that are not available in standard organism annotation databases.

The script is designed for command-line use by lab staff and supports:

- custom GO annotation tables
- gene-level or transcript-level analysis
- input from a differential expression results table or a simple ID list
- `clusterProfiler`, `goseq`, or both
- background universes from the annotation, the DE table, an explicit background file, or an expression matrix
- PDF plots
- robust logging and warning capture
- full and significant-only outputs written into structured folders

## Script

`barley_go_enrichment_custom_v4.R`

## Input files expected

### 1. Annotation file
A tab-delimited file containing at least these columns:

- `BaRTv2 transcript`
- `BaRTv2 gene`
- `GO IDs`
- `GO terms`

Example:

```text
BaRTv2 transcript	BaRTv2 gene	GO IDs	GO terms
BaRT2v18chr1HG000020.1	BaRT2v18chr1HG000020	...	...
```

### 2. Differential expression results file
A file containing significant or full DE results. The script looks for sensible column names such as:

- `target` or `id`
- `contrast`, `comparison`, or similar
- adjusted P-value columns such as `adj.pval`, `padj`, `FDR`, or similar
- log2 fold-change columns such as `log2FC`, `logFC`, or similar
- optional direction labels such as `up.down`

Example:

```text
target,contrast,adj.pval,log2FC,up.down
BaRT2v18chr2HG106390,Live.Barke-HI.Barke,0.0022,-1.6185,down-regulated
```

### 3. Optional ID list file
Instead of a DE file, you can provide a list of IDs.

Supported formats:

- one-column list of IDs
- a file with an ID column plus a set name column, so multiple lists can be analysed in one run

### 4. Optional expression matrix for an expressed background
This can be used to define the background universe from expressed genes or transcripts.

Typical examples:

- `Gene_TPM.csv`
- `Transcript_TPM.csv`

The first column should contain IDs. Sample columns should contain numeric TPM values.

### 5. Optional FASTA file for `goseq`
If `goseq` is run, a FASTA file can be supplied for length-bias correction.

Recommended:

- transcript FASTA for transcript-level analysis
- transcript FASTA is also accepted for gene-level analysis, where transcript lengths are mapped back to genes

## Main features

### Supported analysis modes

- `--id_type gene`
- `--id_type transcript`

### Supported methods

- `--methods clusterprofiler`
- `--methods goseq`
- `--methods both`

### Supported background sources

- `annotation`
- `de_file`
- `background_file`
- `expression_file`

## Significant DE filtering

When a DE file is used, input IDs are selected using:

- `--padj_cutoff`
- `--abs_log2fc_cutoff`

If `--split_direction TRUE`, the script creates up to three sets for each contrast:

- `all`
- `up`
- `down`

## GO significance filtering

The script always writes the full GO result table.

It also writes a second, significant-only version using:

- `--fdr_cutoff`

For `clusterProfiler`, the significant-only file is filtered on:

- `p.adjust <= --fdr_cutoff`

For `goseq`, the significant-only file is filtered on:

- `over_represented_padj <= --fdr_cutoff`
- or `under_represented_padj <= --fdr_cutoff`

## Output structure

Results are written into tool-specific folders inside `--out_dir`.

Example:

```text
<out_dir>/
├── clusterprofiler_GO_enrichment_results/
│   ├── upregulated/
│   │   ├── *.tsv
│   │   ├── *.pdf
│   │   └── significant/
│   │       ├── *_significant_go_enrichment.tsv
│   │       └── *_significant*.pdf
│   ├── downregulated/
│   │   ├── *.tsv
│   │   ├── *.pdf
│   │   └── significant/
│   │       ├── *_significant_go_enrichment.tsv
│   │       └── *_significant*.pdf
│   └── all_significant/
│       ├── *.tsv
│       ├── *.pdf
│       └── significant/
│           ├── *_significant_go_enrichment.tsv
│           └── *_significant*.pdf
├── goseq_GO_enrichment_results/
│   └── ...
├── go_enrichment_run_summary.tsv
├── run_arguments.tsv
├── run_log.txt
├── warnings_log.txt
└── session_info.txt
```

### Notes on the folder names

- `upregulated` contains full and significant-only results for up-regulated IDs
- `downregulated` contains full and significant-only results for down-regulated IDs
- `all_significant` currently contains the combined `all` set results plus its `significant/` subfolder

If you prefer, the `all_significant` folder name can be renamed in a later revision.

## Output files explained

### Full result files

For each analysed contrast or list, the script writes:

- `*_go_enrichment.tsv`
- `*_dotplot.pdf`
- `*_barplot.pdf`

These files contain all returned GO terms.

### Significant-only files

Inside each `significant/` folder, the script writes:

- `*_significant_go_enrichment.tsv`
- significant-only PDF plots when there are any significant terms

### Run summary and logs

- `go_enrichment_run_summary.tsv` summarises every attempted run
- `run_arguments.tsv` records the command-line settings used
- `run_log.txt` records the main steps completed
- `warnings_log.txt` records warnings and recoverable failures
- `session_info.txt` records package and R session information

## Example commands

### 1. Gene-level `clusterProfiler` using a precomputed background file

```bash
Rscript barley_go_enrichment_custom_v4.R \
  --de_file DE_results/sig_DE_genes_list_stats.csv \
  --annotation_file data/BaRT2v18_annotation_GO_only.txt \
  --out_dir go_gene_clusterprofiler \
  --id_type gene \
  --methods clusterprofiler \
  --background_source background_file \
  --background_file data/gene_background_expressed_tpm2_samples2.tsv \
  --split_direction TRUE \
  --make_plots TRUE \
  --fdr_cutoff 0.05
```

### 2. Transcript-level `clusterProfiler` using an expression matrix background

```bash
Rscript barley_go_enrichment_custom_v4.R \
  --de_file DE_results/sig_DE_transcriipts_lists_stats.csv \
  --annotation_file data/BaRT2v18_annotation_GO_only.txt \
  --out_dir go_transcript_clusterprofiler \
  --id_type transcript \
  --methods clusterprofiler \
  --background_source expression_file \
  --expression_file data/Transcript_TPM.csv \
  --expression_min_tpm 2 \
  --expression_min_samples 2 \
  --split_direction TRUE \
  --make_plots TRUE \
  --fdr_cutoff 0.05
```

### 3. Run both `clusterProfiler` and `goseq`

```bash
Rscript barley_go_enrichment_custom_v4.R \
  --de_file DE_results/sig_DE_genes_list_stats.csv \
  --annotation_file data/BaRT2v18_annotation_GO_only.txt \
  --out_dir go_gene_both_methods \
  --id_type gene \
  --methods both \
  --background_source background_file \
  --background_file data/gene_background_expressed_tpm2_samples2.tsv \
  --fasta_file data/BaRT2v18.fa.gz \
  --split_direction TRUE \
  --make_plots TRUE \
  --fdr_cutoff 0.05
```

### 4. Run from a simple list of IDs

```bash
Rscript barley_go_enrichment_custom_v4.R \
  --id_list_file data/my_gene_list.tsv \
  --annotation_file data/BaRT2v18_annotation_GO_only.txt \
  --out_dir go_from_list \
  --id_type gene \
  --methods clusterprofiler \
  --background_source background_file \
  --background_file data/gene_background_expressed_tpm2_samples2.tsv \
  --make_plots TRUE \
  --fdr_cutoff 0.05
```

## Building a background file from a TPM matrix

If you prefer, you can create the expressed background once and then reuse it.

### Gene-level example

```r
expr <- read.csv("data/Gene_TPM.csv", check.names = FALSE)
ids <- trimws(as.character(expr[[1]]))
expr_mat <- as.data.frame(lapply(expr[, -1, drop = FALSE], as.numeric))
keep <- rowSums(expr_mat >= 2, na.rm = TRUE) >= 2
bg_ids <- unique(ids[keep])

write.table(
  x = data.frame(target = bg_ids),
  file = "data/gene_background_expressed_tpm2_samples2.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
```

### Transcript-level example

```r
expr <- read.csv("data/Transcript_TPM.csv", check.names = FALSE)
ids <- trimws(as.character(expr[[1]]))
expr_mat <- as.data.frame(lapply(expr[, -1, drop = FALSE], as.numeric))
keep <- rowSums(expr_mat >= 2, na.rm = TRUE) >= 2
bg_ids <- unique(ids[keep])

write.table(
  x = data.frame(target = bg_ids),
  file = "data/transcript_background_expressed_tpm2_samples2.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
```

## Package requirements

At minimum, the script requires:

- `optparse`
- `clusterProfiler` for `clusterProfiler` analysis
- `goseq` for `goseq` analysis
- plotting packages used by the chosen method

A practical approach is:

- use `--methods clusterprofiler` if `goseq` is difficult to install
- use `--methods both` only when both package stacks are available

## Robustness and failure handling

The script is designed to continue where possible.

- analysis steps are wrapped in `tryCatch()`
- plotting steps are wrapped in `tryCatch()`
- warnings are logged
- failures in one contrast, one direction, or one plotting step should not stop the whole run

This makes it suitable for batch use across many contrasts.

## Troubleshooting

### 1. No background IDs remained after intersecting with the annotation file

This usually means the background IDs do not match the annotation IDs.

Check:

- gene vs transcript mode
- first column of the expression matrix
- whether IDs are in the same BaRT format
- whether a CSV has a blank first column header

A robust workaround is to create an explicit background file first and run with:

- `--background_source background_file`

### 2. `goseq` installation problems

If `goseq` is difficult to install, run only:

```bash
--methods clusterprofiler
```

### 3. No significant GO terms in the `significant/` folder

This means the full enrichment ran, but no terms met:

```text
p.adjust <= --fdr_cutoff
```

or the corresponding `goseq` adjusted P-value threshold.

## Suggested workflow

For most lab runs:

1. create a clean expressed background file once
2. run gene-level enrichment
3. run transcript-level enrichment if needed
4. inspect `go_enrichment_run_summary.tsv` first
5. focus on the `significant/` folders for interpretation and presentation

## Version notes

This README describes the behaviour of `barley_go_enrichment_custom_v4.R`.
