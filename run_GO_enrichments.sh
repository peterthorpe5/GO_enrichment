#!/usr/bin/env bash
#$ -j y
#$ -cwd
#$ -V
#$ -pe smp 2
#$ -jc long
#$ -N GO_enrichment

cd /cluster/db/pthorpe001/Barley_GO_enrichment


conda activate Go_analysis2


IF_NEEDED = """expr <- read.csv("data/Gene_TPM.csv", check.names = FALSE)
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
"""

# with the above made background:

Rscript  GO_enrichment/barley_go_enrichment_custom_v6.R  \
  --de_file DE_results/sig_DE_genes_list_stats.csv \
  --annotation_file data/BaRT2v18_annotation_GO_only.txt \
  --out_dir go_gene_clusterprofiler \
  --id_type gene \
  --methods clusterprofiler \
  --background_source background_file \
  --background_file data/gene_background_expressed_tpm2_samples2.tsv \
  --split_direction TRUE \
  --make_plots TRUE    --fdr_cutoff 0.05

  
# single method_performance

Rscript  GO_enrichment/barley_go_enrichment_custom_v6.R  \
  --de_file DE_results/sig_DE_genes_list_stats.csv \
  --annotation_file data/BaRT2v18_annotation_GO_only.txt \
  --out_dir go_gene_clusterprofiler \
  --id_type gene \
  --methods clusterprofiler \
  --background_source expression_file \
  --expression_file data/Gene_TPM.csv \
  --expression_min_tpm 2 \
  --expression_min_samples 2 \
  --split_direction TRUE \
  --make_plots TRUE --fdr_cutoff 0.05


  Rscript  GO_enrichment/barley_go_enrichment_custom_v6.R  \
  --de_file DE_results/sig_DE_transcripts_lists_stats.csv \
  --annotation_file data/BaRT2v18_annotation_GO_only.txt \
  --out_dir go_transcript_clusterprofiler \
  --id_type transcript \
  --methods clusterprofiler \
  --background_source expression_file \
  --expression_file data/Transcript_TPM.csv \
  --expression_min_tpm 2 \
  --expression_min_samples 2 \
  --split_direction TRUE \
  --make_plots TRUE


  ###### - we need gene lengths for the Goseq

Rscript  GO_enrichment/barley_go_enrichment_custom_v6.R  \
  --de_file data/sig_DE_genes_list_stats.csv \
  --annotation_file data/BaRT2v18_annotation_GO_only.txt \
  --out_dir go_gene_both_methods \
  --id_type gene \
  --methods both \
  --background_source expression_file \
  --expression_file data/Gene_TPM.csv \
  --expression_min_tpm 2 \
  --expression_min_samples 2 \
  --split_direction TRUE \
  --make_plots TRUE    --fdr_cutoff 0.05\
  --fasta_file data/BaRT2v18.fa


Rscript  GO_enrichment/barley_go_enrichment_custom_v6.R  \
  --de_file data/sig_DE_transcriipts_lists_stats.csv \
  --annotation_file data/BaRT2v18_annotation_GO_only.txt \
  --out_dir go_transcript_both_methods \
  --id_type transcript \
  --methods both \
  --background_source expression_file \
  --expression_file data/Transcript_TPM.csv \
  --expression_min_tpm 2 \
  --expression_min_samples 2 \
  --split_direction TRUE \
  --make_plots TRUE    --fdr_cutoff 0.05\
  --fasta_file data/BaRT2v18.fa



