# Custom GO enrichment for barley or any organism with GO annotations

This repository contains an R command-line script for Gene Ontology (GO) enrichment using a custom annotation table. It was developed for barley BaRT-style gene and transcript identifiers, but the workflow is not barley-specific. It can be applied to any organism if you provide an annotation table that links your gene or transcript IDs to GO IDs.

The script is useful when your organism, annotation version, transcriptome assembly, or custom gene set is not available through a standard organism database package.

The script supports:

- custom GO annotation tables
- gene-level or transcript-level analysis
- input from a differential expression results table or a simple ID list
- `clusterProfiler`, `goseq`, or both
- background universes from the annotation, the DE table, an explicit background file, or an expression matrix
- FDR-filtered significant outputs
- GO namespace annotation: `BP`, `MF`, and `CC`
- GO definitions, synonyms, immediate parent terms, and significant ancestor/descendant labels
- zoomed-in most-specific significant GO terms
- zoomed-out broad-parent significant GO terms
- grouping of significant GO terms driven by exactly the same genes
- PDF plots for the full and significant result tables
- tab-separated output tables
- run logs, warning logs, arguments, and session information

## Script

The examples below use:

```text
barley_go_enrichment_custom_v6.R
```

If you have renamed the latest script for compatibility, use your local script name in the commands. The README describes the latest functionality added through the GO-context, broad-parent, most-specific, and gene-set grouping updates.

## What the script does

For each input gene or transcript set, the script:

1. reads the custom annotation table
2. builds a GO-to-gene or GO-to-transcript mapping
3. defines the background universe
4. runs the selected enrichment method or methods
5. writes the full enrichment table
6. writes an FDR-filtered significant table
7. adds GO namespace and GO hierarchy context
8. writes additional interpretation tables for zoomed-in and zoomed-out views
9. optionally groups significant terms driven by the same input genes
10. writes plots, logs, and a run summary

The statistical enrichment results are not changed by the interpretation layers. The extra columns and extra tables are intended to help interpretation.

## Main use cases

### Barley BaRT analysis

Use this with BaRT-style IDs such as:

```text
BaRT2v18chr1HG000020
BaRT2v18chr1HG000020.1
```

### Any organism or custom annotation

You can also use the same script for any organism, for example:

```text
Gene00001
Gene00001.t1
PlasmoGene_0001
Contig123_g00045
TRINITY_DN1000_c0_g1_i1
```

The important requirement is that the same IDs are used consistently across:

- the annotation file
- the DE table or ID list
- the background file or expression matrix
- the FASTA file if using `goseq`

## Input files

All output files are tab-separated. Tab-separated input files are recommended. The script can sometimes fall back to reading comma-separated input if a tab-delimited read produces only one column, but TSV input is safer and easier to debug.

### 1. Annotation file

This file links your features to GO IDs.

The annotation file should contain at least these columns:

```text
transcript
gene
GO IDs
GO terms
```

For BaRT-style annotation files, these column names are also recognised:

```text
BaRTv2 transcript
BaRTv2 gene
GO IDs
GO terms
```

Column matching is case-insensitive and ignores spaces and punctuation. For example, `GO IDs`, `GO_IDs`, and `go ids` are treated similarly.

The script currently looks for the following column-name patterns:

| Role | Recognised names |
|---|---|
| Gene ID | `gene`, `BaRTv2 gene`, `BaRT2 gene` |
| Transcript ID | `transcript`, `BaRTv2 transcript`, `BaRT2 transcript` |
| GO ID list | `GO IDs`, `GO ID` |
| GO term list | `GO terms`, `GO term` |

The `GO IDs` and `GO terms` columns can contain multiple values separated by semicolons or commas.

#### Generic annotation example

```text
gene	transcript	GO IDs	GO terms
Gene00001	Gene00001.t1	GO:0006952;GO:0009611	defense response;response to wounding
Gene00002	Gene00002.t1	GO:0004672;GO:0005524	protein kinase activity;ATP binding
Gene00003	Gene00003.t1	GO:0005737	cytoplasm
```

`GO terms` are useful for readability, but if `GO.db` and `AnnotationDbi` are installed, the script can also retrieve GO term names and definitions from `GO.db`.

#### Barley annotation example

```text
BaRTv2 transcript	BaRTv2 gene	GO IDs	GO terms
BaRT2v18chr1HG000020.1	BaRT2v18chr1HG000020	GO:0005524;GO:0004672	ATP binding;protein kinase activity
BaRT2v18chr1HG000040.1	BaRT2v18chr1HG000040	GO:0005737	cytoplasm
```

#### Annotation file for `goseq` without a FASTA file

If you run `goseq` and do not provide `--fasta_file`, the script tries to calculate feature lengths from annotation coordinates. In that case, the annotation file also needs:

```text
start
end
```

Example:

```text
gene	transcript	GO IDs	GO terms	start	end
Gene00001	Gene00001.t1	GO:0006952	defense response	100	2400
Gene00002	Gene00002.t1	GO:0004672	protein kinase activity	3000	5200
```

For gene-level `goseq`, transcript lengths are collapsed to the maximum transcript length per gene.

### 2. Differential expression results file

Use this when you want the script to select significant IDs from a DE table.

Recommended columns:

```text
target
contrast
adj.pval
log2FC
up.down
```

The script recognises:

| Role | Recognised names |
|---|---|
| Feature ID | `target`, `gene`, `transcript`, `id` |
| Contrast | `contrast` |
| Adjusted P value | `adj.pval` |
| Log2 fold-change | `log2FC` |
| Direction | `up.down` |

The name matching ignores punctuation and case, so `adj.pval` is recognised internally as `adjpval`, and `up.down` as `updown`.

Important: if your DE software writes adjusted P values as `padj`, `FDR`, or another name, rename that column to `adj.pval` before running the script. Otherwise the script may not apply the intended adjusted P-value filter.

#### Generic DE table example

```text
target	contrast	adj.pval	log2FC	up.down
Gene00001	Treatment_vs_Control	0.001	2.4	up-regulated
Gene00002	Treatment_vs_Control	0.020	-1.3	down-regulated
Gene00003	Treatment_vs_Control	0.400	0.2	not-significant
```

The DE table is filtered using:

```text
--padj_cutoff
--abs_log2fc_cutoff
```

If `--split_direction TRUE`, the script creates up to three sets per contrast:

```text
all
up
down
```

If `up.down` is not available, direction is inferred from the sign of `log2FC`.

### 3. Simple ID list file

Use this when you already have a list of genes or transcripts to test.

#### One-list example

```text
target
Gene00001
Gene00002
Gene00003
```

#### Multiple-list example

```text
target	set
Gene00001	root_response
Gene00002	root_response
Gene00003	leaf_response
Gene00004	leaf_response
```

Recognised ID columns are:

```text
target
id
```

Recognised set/list columns are:

```text
set
list
contrast
```

If no set/list column is present, all IDs are analysed as one list named `input_list`.

### 4. Background file

Use this when you want to provide an explicit background universe.

Example:

```text
target
Gene00001
Gene00002
Gene00003
Gene00004
Gene00005
```

The background IDs must be the same type as the analysis mode:

- use gene IDs with `--id_type gene`
- use transcript IDs with `--id_type transcript`

### 5. Expression matrix for an expressed background

Use this when you want the background universe to be all expressed genes or transcripts.

The first column should contain IDs, unless the file has a recognised `target` or `id` column. All other columns should be numeric expression values, usually TPM.

Example:

```text
target	Sample_1	Sample_2	Sample_3
Gene00001	10.2	8.4	0.0
Gene00002	0.1	0.0	0.0
Gene00003	3.2	2.8	5.1
```

The expressed background is defined using:

```text
--expression_min_tpm
--expression_min_samples
```

For example:

```text
--expression_min_tpm 2
--expression_min_samples 2
```

means that a feature must have TPM >= 2 in at least 2 samples to enter the background universe.

### 6. FASTA file for `goseq`

A FASTA file can be supplied for length-bias correction with `goseq`.

Recommended:

- transcript FASTA for transcript-level analysis
- transcript FASTA for gene-level analysis, where transcript lengths are mapped back to genes

The FASTA headers should match the transcript IDs in the annotation file. Anything after the first whitespace in a FASTA header is ignored.

Example:

```text
>Gene00001.t1
ATGGCT...
>Gene00002.t1
ATGACT...
```

If no FASTA is supplied, the script falls back to `start` and `end` columns in the annotation file.

## Main command-line arguments

| Argument | Meaning |
|---|---|
| `--de_file` | DE results file. Optional if `--id_list_file` is supplied. |
| `--id_list_file` | Simple ID-list file. Optional if `--de_file` is supplied. |
| `--annotation_file` | Custom annotation table linking genes/transcripts to GO IDs. |
| `--out_dir` | Output directory. |
| `--id_type` | `gene` or `transcript`. |
| `--methods` | `clusterprofiler`, `goseq`, or `both`. |
| `--background_source` | `annotation`, `de_file`, `background_file`, or `expression_file`. |
| `--background_file` | Explicit background file. |
| `--expression_file` | Expression matrix used to define an expressed background. |
| `--expression_min_tpm` | Minimum TPM threshold for expression background. |
| `--expression_min_samples` | Minimum number of samples meeting TPM threshold. |
| `--fasta_file` | FASTA file for `goseq` length-bias correction. |
| `--padj_cutoff` | Adjusted P-value cutoff used to select significant input IDs from the DE table. |
| `--abs_log2fc_cutoff` | Absolute log2 fold-change cutoff used to select input IDs from the DE table. |
| `--fdr_cutoff` | Adjusted P-value cutoff used to write significant GO result tables. |
| `--min_gs_size` | Minimum GO gene-set size for `clusterProfiler`. |
| `--max_gs_size` | Maximum GO gene-set size for `clusterProfiler`. |
| `--ontology` | `all`, `BP`, `MF`, or `CC`. |
| `--split_direction` | Whether to analyse `all`, `up`, and `down` sets. |
| `--make_plots` | Whether to generate PDF plots. |
| `--top_n_plot` | Number of GO terms to show in plots. |
| `--write_specific_terms` | Write the zoomed-in most-specific significant table. |
| `--write_broad_terms` | Write the zoomed-out broad-parent significant table. |
| `--write_gene_term_groups` | Write the significant gene-set term-group table. |

## Enrichment methods

### `clusterProfiler`

The `clusterProfiler` mode uses `clusterProfiler::enricher()` with the custom `TERM2GENE` and `TERM2NAME` mappings built from the annotation file.

The full result table is filtered only when writing the significant-only files:

```text
p.adjust <= --fdr_cutoff
```

### `goseq`

The `goseq` mode performs GO enrichment with length-bias correction. It needs feature lengths from either:

1. `--fasta_file`, or
2. `start` and `end` columns in the annotation file.

The significant-only `goseq` table is filtered using:

```text
over_represented_padj <= --fdr_cutoff
```

or:

```text
under_represented_padj <= --fdr_cutoff
```

The latest version also adds `geneID` and `Count` columns to `goseq` output where possible, so the gene-set grouping output can be generated for `goseq` as well as `clusterProfiler`.

## GO context columns added to outputs

The output tables include the standard enrichment columns plus additional interpretation columns.

### Namespace columns

```text
ontology
ontology_label
```

These identify whether each term is:

| Code | Meaning |
|---|---|
| `BP` | biological process |
| `MF` | molecular function |
| `CC` | cellular component |

### GO definition and parent-term columns

```text
go_definition
go_synonyms
go_secondary_ids
go_db_ontology
go_immediate_parent_ids
go_immediate_parent_terms
go_immediate_parent_relations
go_immediate_parent_summary
```

These help interpret terms without repeatedly looking them up manually.

For example, a molecular-function esterase term may have an immediate parent such as a broader hydrolase activity. This can help show that several hormone-related terms may actually reflect one broader enzyme-function signal.

### Significant ancestor and descendant columns

```text
go_ancestor_count
significant_descendant_count
has_significant_descendant
most_specific_significant
significant_ancestor_count
has_significant_ancestor
broad_parent_significant
significant_ancestor_terms
significant_descendant_terms
```

These columns are calculated within each significant result set.

They help identify whether a significant GO term is:

- a broad significant parent of other significant terms
- a more specific significant child/descendant term
- a standalone significant term without significant ancestors or descendants in the same result set

## Zoomed-in and zoomed-out result tables

Inside each `significant/` folder, the script writes additional interpretation tables.

### Zoomed-in table: most-specific significant terms

```text
*_significant_most_specific_go_enrichment.tsv
```

This table keeps significant GO terms that do not have significant descendant GO terms in the same result set.

Use this table when you want the more detailed biological interpretation.

### Zoomed-out table: broad-parent significant terms

```text
*_significant_broad_parent_go_enrichment.tsv
```

This table keeps significant GO terms that do not have significant ancestor GO terms in the same result set.

Use this table when you want broader biological themes for summaries, figure legends, or presentations.

### Important interpretation note

GO namespaces are separate:

```text
BP = biological process
MF = molecular function
CC = cellular component
```

A molecular-function term such as an esterase activity will not necessarily sit directly underneath a biological-process term such as a hormone metabolic process. They may be biologically related, but they are not usually parent and child terms in the same GO hierarchy.

The hierarchy columns are therefore best used to understand parent/child relationships within GO, while the gene-set grouping table is best used to detect cases where the same genes are producing several biologically related but different GO terms.

## Gene-set term-group output

The script can write:

```text
*_significant_gene_set_term_groups.tsv
```

This table groups significant GO terms that are driven by exactly the same set of input genes.

This is useful when several significant terms have the same `geneID` values. For example, a small set of genes may be annotated to:

```text
salicylic acid metabolic process
methyl salicylate esterase activity
methyl jasmonate esterase activity
```

The terms have different biological meanings, but the enrichment signal is driven by the same genes. The gene-set term-group output makes this explicit.

Main columns include:

```text
gene_set_group_id
exact_gene_set_shared_by_terms
term_count
gene_set_size
geneID
representative_go_id
representative_description
representative_definition
best_pvalue
best_adjusted_pvalue
ontologies
ontology_labels
go_ids
go_terms
go_definitions
go_synonyms
immediate_parent_terms
immediate_parent_summaries
significant_ancestor_terms
significant_descendant_terms
bp_terms
mf_terms
cc_terms
most_specific_terms
broad_parent_terms
```

To turn this output off:

```bash
--write_gene_term_groups FALSE
```

## Output structure

Results are written into method-specific folders inside `--out_dir`.

Example:

```text
<out_dir>/
├── clusterprofiler_GO_enrichment_results/
│   ├── upregulated/
│   │   ├── *_go_enrichment.tsv
│   │   ├── *_dotplot.pdf
│   │   ├── *_barplot.pdf
│   │   └── significant/
│   │       ├── *_significant_go_enrichment.tsv
│   │       ├── *_significant_most_specific_go_enrichment.tsv
│   │       ├── *_significant_broad_parent_go_enrichment.tsv
│   │       ├── *_significant_gene_set_term_groups.tsv
│   │       ├── *_significant_dotplot.pdf
│   │       └── *_significant_barplot.pdf
│   ├── downregulated/
│   │   └── ...
│   └── all_significant/
│       └── ...
├── goseq_GO_enrichment_results/
│   └── ...
├── go_enrichment_run_summary.tsv
├── run_arguments.tsv
├── run_log.txt
├── warnings_log.txt
└── session_info.txt
```

### Folder names

| Folder | Meaning |
|---|---|
| `upregulated` | Results for up-regulated IDs. |
| `downregulated` | Results for down-regulated IDs. |
| `all_significant` | Results for all selected significant IDs combined. |
| `significant` | FDR-filtered result tables and significant-only plots. |

## Output files explained

### Full result files

For each analysed contrast or list, the script writes:

```text
*_go_enrichment.tsv
*_dotplot.pdf
*_barplot.pdf
```

The full table contains all returned GO terms.

### Significant-only files

Inside each `significant/` folder, the script writes:

```text
*_significant_go_enrichment.tsv
```

This table contains terms passing `--fdr_cutoff`.

### Most-specific significant files

```text
*_significant_most_specific_go_enrichment.tsv
```

This is the zoomed-in table.

### Broad-parent significant files

```text
*_significant_broad_parent_go_enrichment.tsv
```

This is the zoomed-out table.

### Gene-set term-group files

```text
*_significant_gene_set_term_groups.tsv
```

This groups significant terms driven by exactly the same genes.

### Run summary and logs

```text
go_enrichment_run_summary.tsv
run_arguments.tsv
run_log.txt
warnings_log.txt
session_info.txt
```

The run summary includes paths to the main output files, plus counts of returned, significant, and gene-set grouped terms.

## Example commands

### 1. Gene-level `clusterProfiler` using an expression matrix background

```bash
Rscript barley_go_enrichment_custom_v6.R \
  --de_file DE_results/sig_DE_genes_list_stats.tsv \
  --annotation_file data/custom_annotation_GO.tsv \
  --out_dir go_gene_clusterprofiler \
  --id_type gene \
  --methods clusterprofiler \
  --background_source expression_file \
  --expression_file data/Gene_TPM.tsv \
  --expression_min_tpm 2 \
  --expression_min_samples 2 \
  --split_direction TRUE \
  --make_plots TRUE \
  --fdr_cutoff 0.05
```

### 2. Transcript-level `clusterProfiler`

```bash
Rscript barley_go_enrichment_custom_v6.R \
  --de_file DE_results/sig_DE_transcripts_list_stats.tsv \
  --annotation_file data/custom_annotation_GO.tsv \
  --out_dir go_transcript_clusterprofiler \
  --id_type transcript \
  --methods clusterprofiler \
  --background_source expression_file \
  --expression_file data/Transcript_TPM.tsv \
  --expression_min_tpm 2 \
  --expression_min_samples 2 \
  --split_direction TRUE \
  --make_plots TRUE \
  --fdr_cutoff 0.05
```

### 3. Run both `clusterProfiler` and `goseq`

```bash
Rscript barley_go_enrichment_custom_v6.R \
  --de_file DE_results/sig_DE_genes_list_stats.tsv \
  --annotation_file data/custom_annotation_GO.tsv \
  --out_dir go_gene_both_methods \
  --id_type gene \
  --methods both \
  --background_source background_file \
  --background_file data/gene_background_expressed_tpm2_samples2.tsv \
  --fasta_file data/transcripts.fa.gz \
  --split_direction TRUE \
  --make_plots TRUE \
  --fdr_cutoff 0.05
```

### 4. Run from a simple list of IDs

```bash
Rscript barley_go_enrichment_custom_v6.R \
  --id_list_file data/my_gene_list.tsv \
  --annotation_file data/custom_annotation_GO.tsv \
  --out_dir go_from_list \
  --id_type gene \
  --methods clusterprofiler \
  --background_source background_file \
  --background_file data/gene_background_expressed_tpm2_samples2.tsv \
  --make_plots TRUE \
  --fdr_cutoff 0.05
```

### 5. Restrict analysis to one GO namespace

```bash
Rscript barley_go_enrichment_custom_v6.R \
  --de_file DE_results/sig_DE_genes_list_stats.tsv \
  --annotation_file data/custom_annotation_GO.tsv \
  --out_dir go_gene_BP_only \
  --id_type gene \
  --methods clusterprofiler \
  --background_source expression_file \
  --expression_file data/Gene_TPM.tsv \
  --ontology BP \
  --split_direction TRUE \
  --make_plots TRUE \
  --fdr_cutoff 0.05
```

Allowed values are:

```text
all
BP
MF
CC
```

## Building a reusable expressed background file

You can create a clean background file once and reuse it.

### Gene-level example in R

```r
expr <- read.delim(
  file = "data/Gene_TPM.tsv",
  check.names = FALSE
)
ids <- trimws(as.character(expr[[1]]))
expr_mat <- as.data.frame(
  lapply(
    X = expr[, -1, drop = FALSE],
    FUN = as.numeric
  )
)
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

### Transcript-level example in R

```r
expr <- read.delim(
  file = "data/Transcript_TPM.tsv",
  check.names = FALSE
)
ids <- trimws(as.character(expr[[1]]))
expr_mat <- as.data.frame(
  lapply(
    X = expr[, -1, drop = FALSE],
    FUN = as.numeric
  )
)
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

Required packages depend on the selected method.

### Core package

```text
optparse
```

### Required for `clusterProfiler`

```text
clusterProfiler
```

### Required for `goseq`

```text
goseq
```

### Recommended for GO context and `goseq` FASTA handling

```text
GO.db
AnnotationDbi
Biostrings
ggplot2
```

A typical R installation approach is:

```r
install.packages(c("optparse", "ggplot2", "BiocManager"))

BiocManager::install(c(
  "clusterProfiler",
  "goseq",
  "GO.db",
  "AnnotationDbi",
  "Biostrings"
))
```

If `goseq` is difficult to install, run only:

```bash
--methods clusterprofiler
```

## Robustness and failure handling

The script is designed to continue where possible.

- analysis steps are wrapped in `tryCatch()`
- plotting steps are wrapped in `tryCatch()`
- warnings are logged
- failures in one contrast, one direction, or one plotting step should not stop the whole run
- empty result tables are written when no terms are returned
- run arguments and session information are written for reproducibility

## Troubleshooting

### No background IDs remained after intersecting with the annotation file

This usually means the IDs do not match between the background and annotation files.

Check:

- `--id_type gene` versus `--id_type transcript`
- the first column of the expression matrix
- whitespace or version suffixes in IDs
- whether the background uses gene IDs but the annotation analysis uses transcript IDs
- whether the annotation file has the expected `gene` and `transcript` columns

A robust workaround is to create an explicit background file and run with:

```bash
--background_source background_file
```

### The DE file was not filtered as expected

Make sure the adjusted P-value column is named:

```text
adj.pval
```

and the fold-change column is named:

```text
log2FC
```

If your columns are named `padj`, `FDR`, `qvalue`, or similar, rename them before running.

### GO context columns are blank or mostly missing

Install or check:

```text
GO.db
AnnotationDbi
```

The enrichment can still run without rich GO context, but definitions, synonyms, parent terms, and hierarchy-aware columns need these packages.

### `goseq` fails because lengths are missing

Provide one of the following:

1. `--fasta_file` with transcript IDs matching the annotation file, or
2. `start` and `end` columns in the annotation file.

### No significant GO terms in the `significant/` folder

This means the full enrichment ran, but no terms met the selected FDR threshold.

For `clusterProfiler`:

```text
p.adjust <= --fdr_cutoff
```

For `goseq`:

```text
over_represented_padj <= --fdr_cutoff
```

or:

```text
under_represented_padj <= --fdr_cutoff
```

### Several significant terms are driven by the same genes

This is common in GO enrichment.

Use:

```text
*_significant_gene_set_term_groups.tsv
```

This table shows when multiple GO terms are driven by exactly the same input genes. It helps avoid over-interpreting several terms as independent signals.

## Suggested interpretation workflow

For most analyses:

1. inspect `go_enrichment_run_summary.tsv`
2. check the full result table if you need everything returned by the enrichment method
3. use `*_significant_go_enrichment.tsv` as the main formal significant result table
4. use `*_significant_most_specific_go_enrichment.tsv` for more detailed biological interpretation
5. use `*_significant_broad_parent_go_enrichment.tsv` for broader summary themes
6. use `*_significant_gene_set_term_groups.tsv` to identify multiple terms driven by the same genes
7. check `go_definition`, `go_immediate_parent_summary`, `significant_ancestor_terms`, and `significant_descendant_terms` before writing biological conclusions

## Example interpretation of same-gene GO terms

A small gene set may return terms such as:

```text
GO:0009696  salicylic acid metabolic process
GO:0080031  methyl salicylate esterase activity
GO:0080032  methyl jasmonate esterase activity
```

These terms should not be treated as identical.

- `salicylic acid metabolic process` is a biological-process term.
- `methyl salicylate esterase activity` is a molecular-function term.
- `methyl jasmonate esterase activity` is a molecular-function term.

If the same genes drive all three terms, a careful interpretation would be that the result suggests enrichment for hormone-related ester hydrolysis or hormone-conjugate turnover/activation, rather than direct proof that one hormone biosynthetic pathway alone is increased.

The GO context and gene-set grouping columns are designed to make these cases easier to spot.

## Version notes

This README describes the latest GO-context version of the custom enrichment script. Relative to the earlier v4 README, the latest version adds:

- `ontology` and `ontology_label` columns for BP/MF/CC
- GO definitions, synonyms, secondary IDs, and parent-term summaries
- significant ancestor and descendant term labels
- most-specific significant output tables
- broad-parent significant output tables
- significant gene-set term-group output tables
- `geneID` and `Count` support in `goseq` outputs where possible
- richer run-summary paths for the extra outputs

