#!/usr/bin/env Rscript

#' Run GO enrichment for custom barley annotations.
#'
#' This script supports custom Gene Ontology enrichment for barley resources
#' that are not present in standard organism databases. It is designed for
#' simple command-line use by lab staff and accepts BaRT-style annotation files
#' together with either:
#'   1. A differential expression results file.
#'   2. A simple input list of IDs.
#'
#' The script can run clusterProfiler, goseq, or both. It also supports an
#' expressed background universe derived from a TPM matrix.
#'
#' Robustness design:
#'   * Analysis and plotting steps are wrapped in tryCatch().
#'   * Failures are logged as warnings.
#'   * The script continues to the next analysis instead of terminating.
#'   * Plotting failures do not prevent tabular results from being written.
#'
#' Output files are written as tab-separated text files and plots are written
#' as PDF files only.
#'
#' Version note:
#'   * This debug build adds detailed background-universe logging.
#'   * Expression backgrounds now use lapply() rather than apply() to
#'     avoid fragile matrix coercion when reading TPM tables.
#'   * GO ontology context columns are added to output tables.
#'   * Optional most-specific and broad-parent significant tables are written
#'     to help interpret cases where broad parent terms and narrower child
#'     terms are both significant.
#'   * Optional gene-set term-group summaries are written for significant
#'     results. These highlight cases where several significant GO terms are
#'     driven by exactly the same input genes.
#'   * GO definitions, synonyms, immediate parent terms, parent relationship
#'     types, and significant ancestor/descendant term labels are added as
#'     interpretation columns. These columns do not change the statistical
#'     enrichment tests.

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option(
    opt_str = c("--de_file"),
    type = "character",
    default = NULL,
    help = paste(
      "Differential expression results file.",
      "Optional if --id_list_file is supplied."
    )
  ),
  make_option(
    opt_str = c("--id_list_file"),
    type = "character",
    default = NULL,
    help = paste(
      "Optional file containing one or more IDs to test.",
      "Can be one-column, or can include target/id and set/list columns."
    )
  ),
  make_option(
    opt_str = c("--annotation_file"),
    type = "character",
    help = paste(
      "Annotation file with BaRTv2 transcript, BaRTv2 gene, GO IDs,",
      "and GO terms columns."
    )
  ),
  make_option(
    opt_str = c("--out_dir"),
    type = "character",
    default = "go_enrichment_output",
    help = "Output directory. [default: %default]"
  ),
  make_option(
    opt_str = c("--id_type"),
    type = "character",
    default = "gene",
    help = "ID type to analyse: gene or transcript. [default: %default]"
  ),
  make_option(
    opt_str = c("--methods"),
    type = "character",
    default = "both",
    help = paste(
      "Methods to run: clusterprofiler, goseq, or both.",
      "You can also provide a comma-separated list.",
      "[default: %default]"
    )
  ),
  make_option(
    opt_str = c("--background_source"),
    type = "character",
    default = "annotation",
    help = paste(
      "Background source: annotation, de_file, background_file, or expression_file.",
      "[default: %default]"
    )
  ),
  make_option(
    opt_str = c("--background_file"),
    type = "character",
    default = NULL,
    help = "Optional file containing the background universe IDs."
  ),
  make_option(
    opt_str = c("--expression_file"),
    type = "character",
    default = NULL,
    help = paste(
      "Optional expression matrix used to define an expressed background universe.",
      "Suitable for the uploaded Gene TPM.csv or Transcript TPM.csv files."
    )
  ),
  make_option(
    opt_str = c("--expression_min_tpm"),
    type = "double",
    default = 1,
    help = "Minimum TPM for a feature to count as expressed. [default: %default]"
  ),
  make_option(
    opt_str = c("--expression_min_samples"),
    type = "integer",
    default = 1,
    help = paste(
      "Minimum number of samples that must meet the TPM threshold.",
      "[default: %default]"
    )
  ),
  make_option(
    opt_str = c("--fasta_file"),
    type = "character",
    default = NULL,
    help = paste(
      "Optional FASTA file for length bias correction in goseq.",
      "Recommended for transcript-level analysis."
    )
  ),
  make_option(
    opt_str = c("--padj_cutoff"),
    type = "double",
    default = 0.05,
    help = "Adjusted P-value cutoff used to define significant IDs. [default: %default]"
  ),
  make_option(
    opt_str = c("--fdr_cutoff"),
    type = "double",
    default = 0.05,
    help = paste(
      "Adjusted P-value cutoff used to write significant-only GO results.",
      "[default: %default]"
    )
  ),
  make_option(
    opt_str = c("--abs_log2fc_cutoff"),
    type = "double",
    default = 0,
    help = "Absolute log2 fold-change cutoff used to define significant IDs. [default: %default]"
  ),
  make_option(
    opt_str = c("--min_gs_size"),
    type = "integer",
    default = 10,
    help = "Minimum GO term size. [default: %default]"
  ),
  make_option(
    opt_str = c("--max_gs_size"),
    type = "integer",
    default = 500,
    help = "Maximum GO term size. [default: %default]"
  ),
  make_option(
    opt_str = c("--ontology"),
    type = "character",
    default = "all",
    help = "Ontology filter: all, BP, CC, or MF. [default: %default]"
  ),
  make_option(
    opt_str = c("--split_direction"),
    type = "character",
    default = "TRUE",
    help = paste(
      "Whether to analyse up and down separately as well as combined.",
      "TRUE or FALSE. [default: %default]"
    )
  ),
  make_option(
    opt_str = c("--make_plots"),
    type = "character",
    default = "TRUE",
    help = "Whether to generate PDF plots. TRUE or FALSE. [default: %default]"
  ),
  make_option(
    opt_str = c("--top_n_plot"),
    type = "integer",
    default = 20,
    help = "Maximum number of GO terms to include in summary plots. [default: %default]"
  ),
  make_option(
    opt_str = c("--write_specific_terms"),
    type = "character",
    default = "TRUE",
    help = paste(
      "Whether to write an additional significant-most-specific table.",
      "This keeps significant terms that do not have significant descendant",
      "GO terms in the same result set. TRUE or FALSE. [default: %default]"
    )
  ),
  make_option(
    opt_str = c("--write_broad_terms"),
    type = "character",
    default = "TRUE",
    help = paste(
      "Whether to write an additional significant-broad-parent table.",
      "This keeps significant terms that do not have significant ancestor",
      "GO terms in the same result set. TRUE or FALSE. [default: %default]"
    )
  ),
  make_option(
    opt_str = c("--write_gene_term_groups"),
    type = "character",
    default = "TRUE",
    help = paste(
      "Whether to write an additional significant gene-set term-group table.",
      "This groups significant GO terms that are driven by exactly the same",
      "input genes. TRUE or FALSE. [default: %default]"
    )
  )
)

parser <- OptionParser(option_list = option_list)
args <- parse_args(object = parser)

required_args <- c("annotation_file")
missing_args <- required_args[!nzchar(unlist(args[required_args]))]
if (length(x = missing_args) > 0) {
  print_help(object = parser)
  stop(
    sprintf(
      "Missing required argument(s): %s",
      paste(missing_args, collapse = ", ")
    ),
    call. = FALSE
  )
}

if (is.null(x = args$de_file) && is.null(x = args$id_list_file)) {
  print_help(object = parser)
  stop(
    "Please supply either --de_file or --id_list_file.",
    call. = FALSE
  )
}

args$id_type <- tolower(x = args$id_type)
args$background_source <- tolower(x = args$background_source)
args$ontology <- toupper(x = args$ontology)
args$split_direction <- tolower(x = args$split_direction) %in% c("true", "t", "1", "yes", "y")
args$make_plots <- tolower(x = args$make_plots) %in% c("true", "t", "1", "yes", "y")
args$write_specific_terms <- tolower(x = args$write_specific_terms) %in% c("true", "t", "1", "yes", "y")
args$write_broad_terms <- tolower(x = args$write_broad_terms) %in% c("true", "t", "1", "yes", "y")
args$write_gene_term_groups <- tolower(x = args$write_gene_term_groups) %in% c("true", "t", "1", "yes", "y")

parse_methods <- function(method_string) {
  method_string <- tolower(x = trimws(x = method_string))

  if (identical(method_string, "both")) {
    return(c("clusterprofiler", "goseq"))
  }

  methods <- unlist(
    x = strsplit(x = method_string, split = ",", fixed = TRUE),
    use.names = FALSE
  )
  methods <- unique(x = trimws(x = methods))
  methods <- methods[nzchar(x = methods)]
  methods
}

args$methods <- parse_methods(method_string = args$methods)

if (!args$id_type %in% c("gene", "transcript")) {
  stop("--id_type must be 'gene' or 'transcript'.", call. = FALSE)
}

if (!all(args$methods %in% c("clusterprofiler", "goseq"))) {
  stop(
    "--methods must contain only 'clusterprofiler', 'goseq', or 'both'.",
    call. = FALSE
  )
}

if (!args$background_source %in% c("annotation", "de_file", "background_file", "expression_file")) {
  stop(
    paste(
      "--background_source must be 'annotation', 'de_file',",
      "'background_file', or 'expression_file'."
    ),
    call. = FALSE
  )
}

if (!args$ontology %in% c("ALL", "BP", "CC", "MF")) {
  stop("--ontology must be one of: all, BP, CC, MF.", call. = FALSE)
}

read_tabular_file <- function(file_path) {
  df <- utils::read.delim(
    file = file_path,
    header = TRUE,
    sep = "\t",
    quote = "",
    fill = TRUE,
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (ncol(x = df) == 1) {
    df <- utils::read.csv(
      file = file_path,
      header = TRUE,
      quote = '"',
      stringsAsFactors = FALSE,
      check.names = FALSE,
      row.names = NULL
    )
  }

  df
}

normalise_column_names <- function(df) {
  original_names <- colnames(x = df)
  simple_names <- tolower(x = gsub(pattern = "[^A-Za-z0-9]+", replacement = "", x = original_names))
  stats::setNames(object = original_names, nm = simple_names)
}

get_column_name <- function(df, candidates) {
  lookup <- normalise_column_names(df = df)
  hits <- lookup[names(x = lookup) %in% candidates]

  if (length(x = hits) == 0) {
    stop(
      sprintf(
        "Could not find any of the expected columns: %s",
        paste(candidates, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  unname(obj = hits[1])
}

split_multi_value <- function(x) {
  if (length(x = x) == 0 || is.na(x = x) || !nzchar(x = trimws(x = x))) {
    return(character())
  }

  pieces <- unlist(x = strsplit(x = x, split = "[;,]", perl = TRUE), use.names = FALSE)
  pieces <- trimws(x = pieces)
  pieces[nzchar(x = pieces)]
}

safe_dir_name <- function(x) {
  gsub(pattern = "[^A-Za-z0-9._-]+", replacement = "_", x = x)
}

check_packages <- function(methods) {
  required_packages <- c("optparse")

  if ("clusterprofiler" %in% methods) {
    required_packages <- c(required_packages, "clusterProfiler")
  }

  if ("goseq" %in% methods) {
    required_packages <- c(required_packages, "goseq")
  }

  optional_packages <- c("Biostrings", "GO.db", "AnnotationDbi", "ggplot2")
  packages_to_check <- unique(c(required_packages, optional_packages))

  installed <- vapply(
    X = packages_to_check,
    FUN = requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(length = 1)
  )

  missing_required <- required_packages[!installed[required_packages]]

  if (length(x = missing_required) > 0) {
    stop(
      paste(
        "Missing required package(s):",
        paste(missing_required, collapse = ", "),
        "\nInstall from CRAN or Bioconductor before running the script."
      ),
      call. = FALSE
    )
  }

  invisible(installed)
}

parse_annotation <- function(annotation_file, id_type, ontology) {
  annot_df <- read_tabular_file(file_path = annotation_file)

  transcript_col <- get_column_name(
    df = annot_df,
    candidates = c("bartv2transcript", "bart2transcript", "transcript")
  )
  gene_col <- get_column_name(
    df = annot_df,
    candidates = c("bartv2gene", "bart2gene", "gene")
  )
  go_id_col <- get_column_name(
    df = annot_df,
    candidates = c("goids", "goid")
  )
  go_term_col <- get_column_name(
    df = annot_df,
    candidates = c("goterms", "goterm")
  )

  id_col <- if (identical(id_type, "gene")) gene_col else transcript_col

  term2gene_list <- vector(mode = "list", length = nrow(x = annot_df))
  term2name_list <- vector(mode = "list", length = nrow(x = annot_df))

  for (i in seq_len(length.out = nrow(x = annot_df))) {
    feature_id <- trimws(x = as.character(x = annot_df[i, id_col]))
    go_ids <- split_multi_value(x = as.character(x = annot_df[i, go_id_col]))
    go_terms <- split_multi_value(x = as.character(x = annot_df[i, go_term_col]))

    if (!nzchar(x = feature_id) || length(x = go_ids) == 0) {
      next
    }

    term2gene_list[[i]] <- data.frame(
      go_id = go_ids,
      feature_id = rep(x = feature_id, times = length(x = go_ids)),
      stringsAsFactors = FALSE
    )

    if (length(x = go_terms) == length(x = go_ids)) {
      term2name_list[[i]] <- data.frame(
        go_id = go_ids,
        go_term = go_terms,
        stringsAsFactors = FALSE
      )
    } else {
      term2name_list[[i]] <- data.frame(
        go_id = go_ids,
        go_term = rep(x = NA_character_, times = length(x = go_ids)),
        stringsAsFactors = FALSE
      )
    }
  }

  term2gene <- do.call(what = rbind, args = term2gene_list)
  term2name <- do.call(what = rbind, args = term2name_list)

  if (is.null(x = term2gene) || nrow(x = term2gene) == 0) {
    stop("No GO mappings were found in the annotation file.", call. = FALSE)
  }

  term2gene <- unique(term2gene)
  term2name <- unique(term2name)

  if (is.null(x = term2name) || nrow(x = term2name) == 0) {
    term2name <- data.frame(
      go_id = unique(x = term2gene$go_id),
      go_term = NA_character_,
      stringsAsFactors = FALSE
    )
  }

  ontology_df <- NULL

  if (requireNamespace("GO.db", quietly = TRUE) &&
      requireNamespace("AnnotationDbi", quietly = TRUE)) {
    ontology_df <- AnnotationDbi::select(
      x = GO.db::GO.db,
      keys = unique(x = term2gene$go_id),
      columns = c("TERM", "ONTOLOGY"),
      keytype = "GOID"
    )
    ontology_df <- ontology_df[!duplicated(x = ontology_df$GOID), , drop = FALSE]
    colnames(x = ontology_df) <- c("go_id", "go_term_db", "ontology")

    term2name <- merge(
      x = term2name,
      y = ontology_df,
      by = "go_id",
      all.x = TRUE,
      sort = FALSE
    )

    missing_name_idx <- is.na(x = term2name$go_term) | !nzchar(x = term2name$go_term)
    term2name$go_term[missing_name_idx] <- term2name$go_term_db[missing_name_idx]
    term2name$go_term_db <- NULL

    if (!identical(ontology, "ALL")) {
      keep_go <- term2name$go_id[term2name$ontology %in% ontology]
      term2gene <- term2gene[term2gene$go_id %in% keep_go, , drop = FALSE]
      term2name <- term2name[term2name$go_id %in% keep_go, , drop = FALSE]
    }
  }

  list(
    annotation_df = annot_df,
    term2gene = term2gene,
    term2name = term2name,
    transcript_col = transcript_col,
    gene_col = gene_col,
    id_col = id_col,
    ontology_df = ontology_df
  )
}

read_id_vector <- function(file_path) {
  id_df <- read_tabular_file(file_path = file_path)
  lookup <- normalise_column_names(df = id_df)

  if ("target" %in% names(x = lookup)) {
    ids <- id_df[[lookup[["target"]]]]
  } else if ("id" %in% names(x = lookup)) {
    ids <- id_df[[lookup[["id"]]]]
  } else {
    ids <- id_df[[1]]
  }

  ids <- unique(x = trimws(x = as.character(x = ids)))
  ids[nzchar(x = ids)]
}

read_expression_background <- function(expression_file, min_tpm, min_samples) {
  expr_df <- read_tabular_file(file_path = expression_file)
  lookup <- normalise_column_names(df = expr_df)

  id_col <- NULL
  ids <- NULL
  sample_cols <- character()

  if ("target" %in% names(x = lookup)) {
    id_col <- lookup[["target"]]
    ids <- trimws(x = as.character(x = expr_df[[id_col]]))
    sample_cols <- setdiff(x = colnames(x = expr_df), y = id_col)
  } else if ("id" %in% names(x = lookup)) {
    id_col <- lookup[["id"]]
    ids <- trimws(x = as.character(x = expr_df[[id_col]]))
    sample_cols <- setdiff(x = colnames(x = expr_df), y = id_col)
  } else {
    id_col <- colnames(x = expr_df)[1]
    ids <- trimws(x = as.character(x = expr_df[[1]]))
    sample_cols <- colnames(x = expr_df)[-1]
  }

  if (all(is.na(x = ids) | !nzchar(x = ids))) {
    row_ids <- rownames(x = expr_df)
    default_rownames <- as.character(x = seq_len(length.out = nrow(x = expr_df)))

    if (!is.null(x = row_ids) && !all(row_ids %in% default_rownames)) {
      ids <- trimws(x = as.character(x = row_ids))
      sample_cols <- colnames(x = expr_df)
      id_col <- "<row.names>"
    }
  }

  if (length(x = sample_cols) == 0) {
    stop("No sample columns were found in --expression_file.", call. = FALSE)
  }

  expr_mat <- as.data.frame(
    x = lapply(
      X = expr_df[, sample_cols, drop = FALSE],
      FUN = function(x) {
        suppressWarnings(expr = as.numeric(x = x))
      }
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  expressed_idx <- rowSums(x = expr_mat >= min_tpm, na.rm = TRUE) >= min_samples
  ids <- unique(x = ids[expressed_idx])
  ids <- ids[!is.na(x = ids)]
  ids <- ids[nzchar(x = ids)]

  message(sprintf("Expression background file: %s", expression_file))
  message(sprintf("Expression ID column: '%s'", id_col))
  message(sprintf("Number of sample columns: %s", length(x = sample_cols)))
  message(sprintf("Expression threshold: TPM >= %s in >= %s sample(s)", min_tpm, min_samples))
  message(sprintf("Expressed IDs before return: %s", length(x = ids)))
  message(sprintf(
    "First 5 expressed IDs: %s",
    paste(utils::head(x = ids, n = 5), collapse = ", ")
  ))

  ids
}


get_background_ids <- function(
  annotation_info,
  de_df,
  background_source,
  background_file,
  expression_file,
  expression_min_tpm,
  expression_min_samples
) {
  target_col <- NULL

  if (!is.null(x = de_df)) {
    target_col <- get_column_name(
      df = de_df,
      candidates = c("target", "gene", "transcript", "id")
    )
  }

  if (identical(background_source, "annotation")) {
    ids <- unique(x = annotation_info$annotation_df[[annotation_info$id_col]])
  } else if (identical(background_source, "de_file")) {
    if (is.null(x = de_df)) {
      stop("--background_source de_file requires --de_file.", call. = FALSE)
    }
    ids <- unique(x = de_df[[target_col]])
  } else if (identical(background_source, "background_file")) {
    if (is.null(x = background_file)) {
      stop(
        "--background_source was set to background_file but --background_file was not supplied.",
        call. = FALSE
      )
    }
    ids <- read_id_vector(file_path = background_file)
  } else {
    if (is.null(x = expression_file)) {
      stop(
        paste(
          "--background_source was set to expression_file but --expression_file",
          "was not supplied."
        ),
        call. = FALSE
      )
    }
    ids <- read_expression_background(
      expression_file = expression_file,
      min_tpm = expression_min_tpm,
      min_samples = expression_min_samples
    )
  }

  ids <- unique(x = trimws(x = as.character(x = ids)))
  ids[nzchar(x = ids)]
}

extract_lengths_from_fasta <- function(fasta_file, annotation_info, id_type) {
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop(
      "Biostrings is required to read --fasta_file. Please install it first.",
      call. = FALSE
    )
  }

  seqs <- Biostrings::readDNAStringSet(filepath = fasta_file)
  if (requireNamespace("BiocGenerics", quietly = TRUE)) {
    seq_lengths <- BiocGenerics::width(x = seqs)
  } else if (requireNamespace("IRanges", quietly = TRUE)) {
    seq_lengths <- IRanges::width(x = seqs)
  } else {
    seq_lengths <- nchar(x = as.character(x = seqs), type = "chars")
  }
  names(x = seq_lengths) <- sub(pattern = "\\s+.*$", replacement = "", x = names(x = seqs))

  if (identical(id_type, "transcript")) {
    return(seq_lengths)
  }

  mapping_df <- unique(
    annotation_info$annotation_df[, c(annotation_info$transcript_col, annotation_info$gene_col), drop = FALSE]
  )
  colnames(x = mapping_df) <- c("transcript_id", "gene_id")
  mapping_df$seq_length <- unname(obj = seq_lengths[mapping_df$transcript_id])

  gene_lengths <- tapply(
    X = mapping_df$seq_length,
    INDEX = mapping_df$gene_id,
    FUN = function(x) {
      max(x = x, na.rm = TRUE)
    }
  )

  gene_lengths[is.infinite(x = gene_lengths)] <- NA_real_
  gene_lengths
}

extract_lengths_from_annotation <- function(annotation_info, id_type) {
  annot_df <- annotation_info$annotation_df
  start_col <- get_column_name(df = annot_df, candidates = c("start"))
  end_col <- get_column_name(df = annot_df, candidates = c("end"))

  feature_lengths <- abs(x = as.numeric(x = annot_df[[end_col]]) - as.numeric(x = annot_df[[start_col]])) + 1

  if (identical(id_type, "transcript")) {
    names(x = feature_lengths) <- annot_df[[annotation_info$transcript_col]]
    feature_lengths <- tapply(
      X = feature_lengths,
      INDEX = names(x = feature_lengths),
      FUN = function(x) {
        max(x = x, na.rm = TRUE)
      }
    )
    return(feature_lengths)
  }

  gene_ids <- annot_df[[annotation_info$gene_col]]
  gene_lengths <- tapply(
    X = feature_lengths,
    INDEX = gene_ids,
    FUN = function(x) {
      max(x = x, na.rm = TRUE)
    }
  )
  gene_lengths
}

get_length_vector <- function(background_ids, fasta_file, annotation_info, id_type) {
  lengths <- NULL

  if (!is.null(x = fasta_file)) {
    message("Reading lengths from FASTA file.")
    lengths <- extract_lengths_from_fasta(
      fasta_file = fasta_file,
      annotation_info = annotation_info,
      id_type = id_type
    )
  } else {
    message("No FASTA supplied. Falling back to annotation-based genomic spans.")
    lengths <- extract_lengths_from_annotation(
      annotation_info = annotation_info,
      id_type = id_type
    )
  }

  out_lengths <- rep(x = NA_real_, times = length(x = background_ids))
  names(x = out_lengths) <- background_ids
  matched <- background_ids %in% names(x = lengths)
  out_lengths[matched] <- unname(obj = lengths[background_ids[matched]])
  out_lengths
}

prepare_de_sets <- function(de_df, padj_cutoff, abs_log2fc_cutoff, split_direction) {
  target_col <- get_column_name(df = de_df, candidates = c("target", "gene", "transcript", "id"))
  contrast_col <- NULL
  adj_col <- NULL
  lfc_col <- NULL
  direction_col <- NULL

  lookup <- normalise_column_names(df = de_df)

  if ("contrast" %in% names(x = lookup)) {
    contrast_col <- lookup[["contrast"]]
  }

  if ("adjpval" %in% names(x = lookup)) {
    adj_col <- lookup[["adjpval"]]
  }

  if ("log2fc" %in% names(x = lookup)) {
    lfc_col <- lookup[["log2fc"]]
  }

  if ("updown" %in% names(x = lookup)) {
    direction_col <- lookup[["updown"]]
  }

  if (is.null(x = contrast_col)) {
    de_df$.__contrast__ <- "all_contrasts"
    contrast_col <- ".__contrast__"
  }

  contrasts <- unique(x = as.character(x = de_df[[contrast_col]]))
  out <- list()

  for (contrast_name in contrasts) {
    sub_df <- de_df[as.character(x = de_df[[contrast_col]]) %in% contrast_name, , drop = FALSE]

    keep <- rep(x = TRUE, times = nrow(x = sub_df))

    if (!is.null(x = adj_col)) {
      keep <- keep & as.numeric(x = sub_df[[adj_col]]) <= padj_cutoff
    }

    if (!is.null(x = lfc_col)) {
      keep <- keep & abs(x = as.numeric(x = sub_df[[lfc_col]])) >= abs_log2fc_cutoff
    }

    sig_df <- sub_df[keep, , drop = FALSE]
    all_ids <- unique(x = trimws(x = as.character(x = sig_df[[target_col]])))
    all_ids <- all_ids[nzchar(x = all_ids)]

    direction_sets <- list(all = all_ids)

    if (split_direction) {
      if (!is.null(x = direction_col)) {
        up_idx <- grepl(pattern = "up", x = sig_df[[direction_col]], ignore.case = TRUE)
        down_idx <- grepl(pattern = "down", x = sig_df[[direction_col]], ignore.case = TRUE)
      } else if (!is.null(x = lfc_col)) {
        up_idx <- as.numeric(x = sig_df[[lfc_col]]) > 0
        down_idx <- as.numeric(x = sig_df[[lfc_col]]) < 0
      } else {
        up_idx <- rep(x = FALSE, times = nrow(x = sig_df))
        down_idx <- rep(x = FALSE, times = nrow(x = sig_df))
      }

      direction_sets$up <- unique(x = trimws(x = as.character(x = sig_df[[target_col]][up_idx])))
      direction_sets$down <- unique(x = trimws(x = as.character(x = sig_df[[target_col]][down_idx])))

      direction_sets$up <- direction_sets$up[nzchar(x = direction_sets$up)]
      direction_sets$down <- direction_sets$down[nzchar(x = direction_sets$down)]
    }

    out[[contrast_name]] <- direction_sets
  }

  out
}

prepare_id_list_sets <- function(id_list_file) {
  list_df <- read_tabular_file(file_path = id_list_file)
  lookup <- normalise_column_names(df = list_df)

  target_col <- NULL
  set_col <- NULL

  if ("target" %in% names(x = lookup)) {
    target_col <- lookup[["target"]]
  } else if ("id" %in% names(x = lookup)) {
    target_col <- lookup[["id"]]
  } else {
    target_col <- colnames(x = list_df)[1]
  }

  if ("set" %in% names(x = lookup)) {
    set_col <- lookup[["set"]]
  } else if ("list" %in% names(x = lookup)) {
    set_col <- lookup[["list"]]
  } else if ("contrast" %in% names(x = lookup)) {
    set_col <- lookup[["contrast"]]
  }

  if (is.null(x = set_col)) {
    ids <- unique(x = trimws(x = as.character(x = list_df[[target_col]])))
    ids <- ids[nzchar(x = ids)]
    return(list(input_list = list(all = ids)))
  }

  set_names <- unique(x = as.character(x = list_df[[set_col]]))
  out <- list()

  for (set_name in set_names) {
    sub_df <- list_df[as.character(x = list_df[[set_col]]) %in% set_name, , drop = FALSE]
    ids <- unique(x = trimws(x = as.character(x = sub_df[[target_col]])))
    ids <- ids[nzchar(x = ids)]
    out[[set_name]] <- list(all = ids)
  }

  out
}

run_clusterprofiler <- function(
  sig_ids,
  background_ids,
  term2gene,
  term2name,
  min_gs_size,
  max_gs_size
) {
  filtered_term2gene <- term2gene[term2gene$feature_id %in% background_ids, , drop = FALSE]
  filtered_term2name <- term2name[term2name$go_id %in% filtered_term2gene$go_id, , drop = FALSE]

  if (length(x = sig_ids) == 0 || nrow(x = filtered_term2gene) == 0) {
    return(list(result_df = data.frame(), result_obj = NULL))
  }

  result_obj <- clusterProfiler::enricher(
    gene = sig_ids,
    universe = background_ids,
    TERM2GENE = filtered_term2gene[, c("go_id", "feature_id")],
    TERM2NAME = filtered_term2name[, c("go_id", "go_term")],
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    minGSSize = min_gs_size,
    maxGSSize = max_gs_size
  )

  if (is.null(x = result_obj)) {
    return(list(result_df = data.frame(), result_obj = NULL))
  }

  result_df <- as.data.frame(x = result_obj)
  list(result_df = result_df, result_obj = result_obj)
}

run_goseq <- function(
  sig_ids,
  background_ids,
  term2gene,
  term2name,
  length_vector
) {
  gene_vector <- as.integer(x = background_ids %in% sig_ids)
  names(x = gene_vector) <- background_ids

  gene2cat <- term2gene[term2gene$feature_id %in% background_ids, c("feature_id", "go_id"), drop = FALSE]
  colnames(x = gene2cat) <- c("gene", "category")

  if (length(x = sig_ids) == 0 || nrow(x = gene2cat) == 0) {
    return(list(result_df = data.frame(), pwf = NULL))
  }

  pwf <- goseq::nullp(
    DEgenes = gene_vector,
    bias.data = length_vector,
    plot.fit = FALSE
  )

  result_df <- goseq::goseq(
    pwf = pwf,
    gene2cat = gene2cat,
    use_genes_without_cat = TRUE,
    method = "Wallenius"
  )

  result_df$over_represented_padj <- stats::p.adjust(
    p = result_df$over_represented_pvalue,
    method = "BH"
  )
  result_df$under_represented_padj <- stats::p.adjust(
    p = result_df$under_represented_pvalue,
    method = "BH"
  )

  result_df <- merge(
    x = result_df,
    y = unique(x = term2name[, intersect(
      x = c("go_id", "go_term", "ontology"),
      y = colnames(x = term2name)
    ), drop = FALSE]),
    by.x = "category",
    by.y = "go_id",
    all.x = TRUE,
    sort = FALSE
  )

  de_gene2cat <- gene2cat[gene2cat$gene %in% sig_ids, , drop = FALSE]
  if (nrow(x = de_gene2cat) > 0) {
    gene_lists <- stats::aggregate(
      x = de_gene2cat$gene,
      by = list(category = de_gene2cat$category),
      FUN = function(x) {
        paste(sort(x = unique(x = as.character(x))), collapse = "/")
      }
    )
    colnames(x = gene_lists) <- c("category", "geneID")
    gene_lists$Count <- vapply(
      X = strsplit(x = gene_lists$geneID, split = "/", fixed = TRUE),
      FUN = length,
      FUN.VALUE = integer(length = 1)
    )

    result_df <- merge(
      x = result_df,
      y = gene_lists,
      by = "category",
      all.x = TRUE,
      sort = FALSE
    )
  } else {
    result_df$geneID <- NA_character_
    result_df$Count <- 0L
  }

  list(result_df = result_df, pwf = pwf)
}


get_ontology_label <- function(ontology_code) {
  labels <- c(
    BP = "biological_process",
    MF = "molecular_function",
    CC = "cellular_component"
  )
  out <- unname(obj = labels[ontology_code])
  out[is.na(x = out)] <- NA_character_
  out
}


get_go_term_names <- function(go_ids) {
  go_ids <- as.character(x = go_ids)
  out <- stats::setNames(
    object = rep(x = NA_character_, times = length(x = go_ids)),
    nm = go_ids
  )

  if (!requireNamespace("GO.db", quietly = TRUE) ||
      !requireNamespace("AnnotationDbi", quietly = TRUE)) {
    return(out)
  }

  for (go_id in unique(x = go_ids)) {
    if (is.na(x = go_id) || !nzchar(x = go_id)) {
      next
    }

    term_object <- tryCatch(
      expr = GO.db::GOTERM[[go_id]],
      error = function(e) NULL
    )

    if (is.null(x = term_object)) {
      next
    }

    term_name <- tryCatch(
      expr = AnnotationDbi::Term(term_object),
      error = function(e) NA_character_
    )
    term_name <- as.character(x = term_name)
    term_name <- term_name[!is.na(x = term_name) & nzchar(x = term_name)]

    if (length(x = term_name) > 0) {
      out[names(x = out) %in% go_id] <- term_name[[1]]
    }
  }

  out
}

format_go_term_labels <- function(go_ids, term_lookup = NULL) {
  go_ids <- unique(x = as.character(x = go_ids))
  go_ids <- go_ids[!is.na(x = go_ids) & nzchar(x = go_ids)]

  if (length(x = go_ids) == 0) {
    return("")
  }

  if (is.null(x = term_lookup)) {
    terms <- get_go_term_names(go_ids = go_ids)
  } else {
    terms <- as.character(x = term_lookup[go_ids])
    names(x = terms) <- go_ids
    missing_terms <- is.na(x = terms) | !nzchar(x = terms)
    if (any(missing_terms)) {
      terms[missing_terms] <- get_go_term_names(go_ids = go_ids[missing_terms])
    }
  }

  terms[is.na(x = terms) | !nzchar(x = terms)] <- "unknown_term"
  paste(paste(go_ids, terms, sep = " = "), collapse = " | ")
}

get_go_basic_info <- function(go_ids) {
  go_ids <- unique(x = as.character(x = go_ids))
  go_ids <- go_ids[!is.na(x = go_ids) & nzchar(x = go_ids)]

  out <- data.frame(
    go_id = go_ids,
    go_definition = NA_character_,
    go_synonyms = NA_character_,
    go_secondary_ids = NA_character_,
    go_db_ontology = NA_character_,
    stringsAsFactors = FALSE
  )

  if (nrow(x = out) == 0 ||
      !requireNamespace("GO.db", quietly = TRUE) ||
      !requireNamespace("AnnotationDbi", quietly = TRUE)) {
    return(out)
  }

  for (i in seq_len(length.out = nrow(x = out))) {
    go_id <- out$go_id[[i]]
    term_object <- tryCatch(
      expr = GO.db::GOTERM[[go_id]],
      error = function(e) NULL
    )

    if (is.null(x = term_object)) {
      next
    }

    out$go_definition[[i]] <- collapse_values(tryCatch(
      expr = AnnotationDbi::Definition(term_object),
      error = function(e) NA_character_
    ))
    out$go_synonyms[[i]] <- collapse_values(tryCatch(
      expr = AnnotationDbi::Synonym(term_object),
      error = function(e) NA_character_
    ))
    out$go_secondary_ids[[i]] <- collapse_values(tryCatch(
      expr = AnnotationDbi::Secondary(term_object),
      error = function(e) NA_character_
    ))
    out$go_db_ontology[[i]] <- collapse_values(tryCatch(
      expr = AnnotationDbi::Ontology(term_object),
      error = function(e) NA_character_
    ))
  }

  out
}

get_go_parents <- function(go_ids, ontology_codes) {
  out <- stats::setNames(
    object = vector(mode = "list", length = length(x = go_ids)),
    nm = go_ids
  )

  for (go_id in names(x = out)) {
    out[[go_id]] <- list(ids = character(), relations = character())
  }

  if (!requireNamespace("GO.db", quietly = TRUE) ||
      !requireNamespace("AnnotationDbi", quietly = TRUE)) {
    return(out)
  }

  env_lookup <- list(
    BP = GO.db::GOBPPARENTS,
    MF = GO.db::GOMFPARENTS,
    CC = GO.db::GOCCPARENTS
  )

  for (ontology_code in names(x = env_lookup)) {
    these_ids <- go_ids[ontology_codes %in% ontology_code]
    these_ids <- these_ids[
      these_ids %in% AnnotationDbi::mappedkeys(x = env_lookup[[ontology_code]])
    ]

    if (length(x = these_ids) == 0) {
      next
    }

    parent_lookup <- as.list(x = env_lookup[[ontology_code]][these_ids])

    for (go_id in names(x = parent_lookup)) {
      value <- parent_lookup[[go_id]]
      if (is.null(x = value) || length(x = value) == 0) {
        next
      }

      raw_parent_ids <- as.character(x = unname(obj = value))
      raw_relations <- as.character(x = names(x = value))
      raw_relations[is.na(x = raw_relations) | !nzchar(x = raw_relations)] <- "unknown_relation"

      keep <- !is.na(x = raw_parent_ids) & nzchar(x = raw_parent_ids)
      raw_parent_ids <- raw_parent_ids[keep]
      raw_relations <- raw_relations[keep]

      if (length(x = raw_parent_ids) == 0) {
        next
      }

      dedup_index <- !duplicated(x = raw_parent_ids)
      out[[go_id]] <- list(
        ids = raw_parent_ids[dedup_index],
        relations = raw_relations[dedup_index]
      )
    }
  }

  out
}

format_go_parent_summary <- function(parent_entry) {
  if (is.null(x = parent_entry) || length(x = parent_entry$ids) == 0) {
    return("")
  }

  parent_ids <- parent_entry$ids
  parent_terms <- get_go_term_names(go_ids = parent_ids)
  parent_terms[is.na(x = parent_terms) | !nzchar(x = parent_terms)] <- "unknown_term"
  relations <- parent_entry$relations

  if (length(x = relations) != length(x = parent_ids)) {
    relations <- rep(x = "unknown_relation", times = length(x = parent_ids))
  }

  paste(
    sprintf(
      "%s [%s] %s",
      parent_ids,
      relations,
      as.character(x = parent_terms[parent_ids])
    ),
    collapse = " | "
  )
}

get_go_ancestors <- function(go_ids, ontology_codes) {
  out <- stats::setNames(
    object = vector(mode = "list", length = length(x = go_ids)),
    nm = go_ids
  )

  if (!requireNamespace("GO.db", quietly = TRUE) ||
      !requireNamespace("AnnotationDbi", quietly = TRUE)) {
    return(out)
  }

  env_lookup <- list(
    BP = GO.db::GOBPANCESTOR,
    MF = GO.db::GOMFANCESTOR,
    CC = GO.db::GOCCANCESTOR
  )

  for (ontology_code in names(x = env_lookup)) {
    these_ids <- go_ids[ontology_codes %in% ontology_code]
    these_ids <- these_ids[these_ids %in% AnnotationDbi::mappedkeys(x = env_lookup[[ontology_code]])]

    if (length(x = these_ids) == 0) {
      next
    }

    ancestor_lookup <- as.list(x = env_lookup[[ontology_code]][these_ids])
    ancestor_list <- stats::setNames(
      object = lapply(
        X = these_ids,
        FUN = function(go_id) {
          value <- ancestor_lookup[[go_id]]
          if (is.null(x = value)) {
            return(NA_character_)
          }
          value
        }
      ),
      nm = these_ids
    )

    for (go_id in names(x = ancestor_list)) {
      ancestors <- unique(x = as.character(x = ancestor_list[[go_id]]))
      ancestors <- ancestors[!is.na(x = ancestors) & nzchar(x = ancestors)]
      out[[go_id]] <- ancestors
    }
  }

  out
}

add_go_context_columns <- function(result_df, method, term2name) {
  if (nrow(x = result_df) == 0) {
    return(result_df)
  }

  id_col <- if (identical(method, "clusterprofiler")) "ID" else "category"

  if (!id_col %in% colnames(x = result_df)) {
    return(result_df)
  }

  result_df$.original_order <- seq_len(length.out = nrow(x = result_df))

  context_cols <- c(
    "ontology",
    "ontology_label",
    "go_definition",
    "go_synonyms",
    "go_secondary_ids",
    "go_db_ontology",
    "go_immediate_parent_ids",
    "go_immediate_parent_terms",
    "go_immediate_parent_relations",
    "go_immediate_parent_summary",
    "go_ancestor_count",
    "significant_descendant_count",
    "has_significant_descendant",
    "most_specific_significant",
    "significant_ancestor_count",
    "has_significant_ancestor",
    "broad_parent_significant",
    "significant_ancestor_terms",
    "significant_descendant_terms"
  )
  result_df <- result_df[, !colnames(x = result_df) %in% context_cols, drop = FALSE]

  meta_cols <- intersect(
    x = c("go_id", "ontology"),
    y = colnames(x = term2name)
  )
  meta_df <- unique(x = term2name[, meta_cols, drop = FALSE])

  result_df <- merge(
    x = result_df,
    y = meta_df,
    by.x = id_col,
    by.y = "go_id",
    all.x = TRUE,
    sort = FALSE
  )
  result_df <- result_df[order(result_df$.original_order), , drop = FALSE]
  result_df$.original_order <- NULL

  if (!"ontology" %in% colnames(x = result_df)) {
    result_df$ontology <- NA_character_
  }

  result_df$ontology_label <- get_ontology_label(ontology_code = result_df$ontology)

  go_info_df <- get_go_basic_info(
    go_ids = as.character(x = result_df[[id_col]])
  )
  if (nrow(x = go_info_df) > 0) {
    result_df$.original_order <- seq_len(length.out = nrow(x = result_df))
    result_df <- merge(
      x = result_df,
      y = go_info_df,
      by.x = id_col,
      by.y = "go_id",
      all.x = TRUE,
      sort = FALSE
    )
    result_df <- result_df[order(result_df$.original_order), , drop = FALSE]
    result_df$.original_order <- NULL
  }

  parent_list <- get_go_parents(
    go_ids = as.character(x = result_df[[id_col]]),
    ontology_codes = as.character(x = result_df$ontology)
  )

  ordered_ids <- as.character(x = result_df[[id_col]])
  result_df$go_immediate_parent_ids <- vapply(
    X = parent_list[ordered_ids],
    FUN = function(parent_entry) {
      paste(parent_entry$ids, collapse = " | ")
    },
    FUN.VALUE = character(length = 1)
  )
  result_df$go_immediate_parent_terms <- vapply(
    X = parent_list[ordered_ids],
    FUN = function(parent_entry) {
      if (length(x = parent_entry$ids) == 0) {
        return("")
      }
      parent_terms <- get_go_term_names(go_ids = parent_entry$ids)
      collapse_values(x = as.character(x = parent_terms[parent_entry$ids]))
    },
    FUN.VALUE = character(length = 1)
  )
  result_df$go_immediate_parent_relations <- vapply(
    X = parent_list[ordered_ids],
    FUN = function(parent_entry) {
      paste(parent_entry$relations, collapse = " | ")
    },
    FUN.VALUE = character(length = 1)
  )
  result_df$go_immediate_parent_summary <- vapply(
    X = parent_list[ordered_ids],
    FUN = format_go_parent_summary,
    FUN.VALUE = character(length = 1)
  )

  ancestor_list <- get_go_ancestors(
    go_ids = as.character(x = result_df[[id_col]]),
    ontology_codes = as.character(x = result_df$ontology)
  )

  result_df$go_ancestor_count <- vapply(
    X = ancestor_list[as.character(x = result_df[[id_col]])],
    FUN = length,
    FUN.VALUE = integer(length = 1)
  )

  result_df$significant_descendant_count <- 0L
  result_df$has_significant_descendant <- FALSE
  result_df$most_specific_significant <- NA
  result_df$significant_ancestor_count <- 0L
  result_df$has_significant_ancestor <- FALSE
  result_df$broad_parent_significant <- NA
  result_df$significant_ancestor_terms <- ""
  result_df$significant_descendant_terms <- ""

  result_df
}

mark_go_hierarchy_significance <- function(significant_df, method) {
  if (nrow(x = significant_df) == 0) {
    return(significant_df)
  }

  id_col <- if (identical(method, "clusterprofiler")) "ID" else "category"
  term_col <- get_result_term_column(method = method)

  if (!id_col %in% colnames(x = significant_df) ||
      !"ontology" %in% colnames(x = significant_df)) {
    return(significant_df)
  }

  go_ids <- as.character(x = significant_df[[id_col]])
  ancestor_list <- get_go_ancestors(
    go_ids = go_ids,
    ontology_codes = as.character(x = significant_df$ontology)
  )

  descendant_count <- integer(length = length(x = go_ids))
  ancestor_count <- integer(length = length(x = go_ids))
  names(x = descendant_count) <- go_ids
  names(x = ancestor_count) <- go_ids

  significant_ancestor_ids <- stats::setNames(
    object = vector(mode = "list", length = length(x = go_ids)),
    nm = go_ids
  )
  significant_descendant_ids <- stats::setNames(
    object = vector(mode = "list", length = length(x = go_ids)),
    nm = go_ids
  )

  for (go_id in go_ids) {
    significant_ancestor_ids[[go_id]] <- character()
    significant_descendant_ids[[go_id]] <- character()
  }

  for (child_id in go_ids) {
    ancestors <- ancestor_list[[child_id]]
    if (length(x = ancestors) == 0) {
      next
    }

    significant_ancestors <- intersect(x = ancestors, y = go_ids)
    if (length(x = significant_ancestors) == 0) {
      next
    }

    ancestor_count[child_id] <- length(x = significant_ancestors)
    significant_ancestor_ids[[child_id]] <- significant_ancestors

    for (ancestor_id in significant_ancestors) {
      descendant_count[ancestor_id] <- descendant_count[ancestor_id] + 1L
      significant_descendant_ids[[ancestor_id]] <- unique(
        x = c(significant_descendant_ids[[ancestor_id]], child_id)
      )
    }
  }

  term_lookup <- NULL
  if (term_col %in% colnames(x = significant_df)) {
    term_lookup <- stats::setNames(
      object = as.character(x = significant_df[[term_col]]),
      nm = go_ids
    )
  }

  significant_df$significant_descendant_count <- as.integer(x = descendant_count[go_ids])
  significant_df$has_significant_descendant <- significant_df$significant_descendant_count > 0L
  significant_df$most_specific_significant <- !significant_df$has_significant_descendant
  significant_df$significant_ancestor_count <- as.integer(x = ancestor_count[go_ids])
  significant_df$has_significant_ancestor <- significant_df$significant_ancestor_count > 0L
  significant_df$broad_parent_significant <- !significant_df$has_significant_ancestor
  significant_df$significant_ancestor_terms <- vapply(
    X = significant_ancestor_ids[go_ids],
    FUN = format_go_term_labels,
    FUN.VALUE = character(length = 1),
    term_lookup = term_lookup
  )
  significant_df$significant_descendant_terms <- vapply(
    X = significant_descendant_ids[go_ids],
    FUN = format_go_term_labels,
    FUN.VALUE = character(length = 1),
    term_lookup = term_lookup
  )

  significant_df
}

mark_most_specific_significant <- function(significant_df, method) {
  mark_go_hierarchy_significance(
    significant_df = significant_df,
    method = method
  )
}

filter_most_specific_significant <- function(significant_df, method) {
  marked_df <- mark_go_hierarchy_significance(
    significant_df = significant_df,
    method = method
  )

  if (!"most_specific_significant" %in% colnames(x = marked_df)) {
    return(marked_df)
  }

  marked_df[marked_df$most_specific_significant %in% TRUE, , drop = FALSE]
}

filter_broad_parent_significant <- function(significant_df, method) {
  marked_df <- mark_go_hierarchy_significance(
    significant_df = significant_df,
    method = method
  )

  if (!"broad_parent_significant" %in% colnames(x = marked_df)) {
    return(marked_df)
  }

  marked_df[marked_df$broad_parent_significant %in% TRUE, , drop = FALSE]
}

get_result_id_column <- function(method) {
  if (identical(method, "clusterprofiler")) {
    return("ID")
  }

  "category"
}

get_result_term_column <- function(method) {
  if (identical(method, "clusterprofiler")) {
    return("Description")
  }

  "go_term"
}

split_gene_id_string <- function(gene_string) {
  if (length(x = gene_string) == 0 ||
      is.na(x = gene_string) ||
      !nzchar(x = trimws(x = gene_string))) {
    return(character())
  }

  genes <- unlist(
    x = strsplit(x = gene_string, split = "/", fixed = TRUE),
    use.names = FALSE
  )
  genes <- trimws(x = genes)
  genes <- genes[nzchar(x = genes)]
  sort(x = unique(x = genes))
}

get_p_adjust_column <- function(result_df, method) {
  if (identical(method, "clusterprofiler") &&
      "p.adjust" %in% colnames(x = result_df)) {
    return("p.adjust")
  }

  if ("over_represented_padj" %in% colnames(x = result_df)) {
    return("over_represented_padj")
  }

  NA_character_
}

get_p_value_column <- function(result_df, method) {
  if (identical(method, "clusterprofiler") &&
      "pvalue" %in% colnames(x = result_df)) {
    return("pvalue")
  }

  if ("over_represented_pvalue" %in% colnames(x = result_df)) {
    return("over_represented_pvalue")
  }

  NA_character_
}

collapse_values <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x = x) & nzchar(x = trimws(x = x))]
  paste(unique(x = x), collapse = " | ")
}

make_empty_gene_term_group_result <- function() {
  data.frame(
    gene_set_group_id = character(),
    exact_gene_set_shared_by_terms = logical(),
    term_count = integer(),
    gene_set_size = integer(),
    geneID = character(),
    representative_go_id = character(),
    representative_description = character(),
    representative_definition = character(),
    best_pvalue = numeric(),
    best_adjusted_pvalue = numeric(),
    ontologies = character(),
    ontology_labels = character(),
    go_ids = character(),
    go_terms = character(),
    go_definitions = character(),
    go_synonyms = character(),
    immediate_parent_terms = character(),
    immediate_parent_summaries = character(),
    significant_ancestor_terms = character(),
    significant_descendant_terms = character(),
    bp_terms = character(),
    mf_terms = character(),
    cc_terms = character(),
    most_specific_terms = character(),
    broad_parent_terms = character(),
    stringsAsFactors = FALSE
  )
}

make_gene_set_term_groups <- function(significant_df, method) {
  if (nrow(x = significant_df) == 0 ||
      !"geneID" %in% colnames(x = significant_df)) {
    return(make_empty_gene_term_group_result())
  }

  id_col <- get_result_id_column(method = method)
  term_col <- get_result_term_column(method = method)

  if (!id_col %in% colnames(x = significant_df) ||
      !term_col %in% colnames(x = significant_df)) {
    return(make_empty_gene_term_group_result())
  }

  significant_df$gene_set_signature <- vapply(
    X = significant_df$geneID,
    FUN = function(gene_string) {
      paste(split_gene_id_string(gene_string = gene_string), collapse = "/")
    },
    FUN.VALUE = character(length = 1)
  )

  significant_df <- significant_df[nzchar(x = significant_df$gene_set_signature), , drop = FALSE]

  if (nrow(x = significant_df) == 0) {
    return(make_empty_gene_term_group_result())
  }

  padj_col <- get_p_adjust_column(result_df = significant_df, method = method)
  pval_col <- get_p_value_column(result_df = significant_df, method = method)

  group_keys <- unique(x = significant_df$gene_set_signature)
  group_rows <- vector(mode = "list", length = length(x = group_keys))

  for (i in seq_along(along.with = group_keys)) {
    group_key <- group_keys[[i]]
    group_df <- significant_df[
      significant_df$gene_set_signature %in% group_key,
      ,
      drop = FALSE
    ]

    if (!is.na(x = padj_col)) {
      best_order <- order(group_df[[padj_col]], na.last = TRUE)
    } else if (!is.na(x = pval_col)) {
      best_order <- order(group_df[[pval_col]], na.last = TRUE)
    } else {
      best_order <- seq_len(length.out = nrow(x = group_df))
    }

    representative <- group_df[best_order[[1]], , drop = FALSE]
    genes <- split_gene_id_string(gene_string = group_key)

    ontology_values <- if ("ontology" %in% colnames(x = group_df)) {
      collapse_values(x = group_df$ontology)
    } else {
      ""
    }

    ontology_label_values <- if ("ontology_label" %in% colnames(x = group_df)) {
      collapse_values(x = group_df$ontology_label)
    } else {
      ""
    }

    best_padjust <- if (!is.na(x = padj_col)) {
      suppressWarnings(min(group_df[[padj_col]], na.rm = TRUE))
    } else {
      NA_real_
    }
    if (is.infinite(x = best_padjust)) {
      best_padjust <- NA_real_
    }

    best_pvalue <- if (!is.na(x = pval_col)) {
      suppressWarnings(min(group_df[[pval_col]], na.rm = TRUE))
    } else {
      NA_real_
    }
    if (is.infinite(x = best_pvalue)) {
      best_pvalue <- NA_real_
    }

    term_labels <- paste(
      as.character(x = group_df[[id_col]]),
      as.character(x = group_df[[term_col]]),
      sep = " = "
    )

    go_definition_values <- if ("go_definition" %in% colnames(x = group_df)) {
      collapse_values(x = paste(
        as.character(x = group_df[[id_col]]),
        as.character(x = group_df$go_definition),
        sep = " = "
      ))
    } else {
      ""
    }

    go_synonym_values <- if ("go_synonyms" %in% colnames(x = group_df)) {
      collapse_values(x = paste(
        as.character(x = group_df[[id_col]]),
        as.character(x = group_df$go_synonyms),
        sep = " = "
      ))
    } else {
      ""
    }

    immediate_parent_values <- if ("go_immediate_parent_terms" %in% colnames(x = group_df)) {
      collapse_values(x = group_df$go_immediate_parent_terms)
    } else {
      ""
    }

    immediate_parent_summary_values <- if ("go_immediate_parent_summary" %in% colnames(x = group_df)) {
      collapse_values(x = group_df$go_immediate_parent_summary)
    } else {
      ""
    }

    significant_ancestor_values <- if ("significant_ancestor_terms" %in% colnames(x = group_df)) {
      collapse_values(x = group_df$significant_ancestor_terms)
    } else {
      ""
    }

    significant_descendant_values <- if ("significant_descendant_terms" %in% colnames(x = group_df)) {
      collapse_values(x = group_df$significant_descendant_terms)
    } else {
      ""
    }

    representative_definition <- if ("go_definition" %in% colnames(x = representative)) {
      as.character(x = representative$go_definition[[1]])
    } else {
      ""
    }

    bp_terms <- character()
    mf_terms <- character()
    cc_terms <- character()
    if ("ontology" %in% colnames(x = group_df)) {
      bp_terms <- term_labels[group_df$ontology %in% "BP"]
      mf_terms <- term_labels[group_df$ontology %in% "MF"]
      cc_terms <- term_labels[group_df$ontology %in% "CC"]
    }

    most_specific_terms <- character()
    if ("most_specific_significant" %in% colnames(x = group_df)) {
      most_specific_terms <- term_labels[group_df$most_specific_significant %in% TRUE]
    }

    broad_parent_terms <- character()
    if ("broad_parent_significant" %in% colnames(x = group_df)) {
      broad_parent_terms <- term_labels[group_df$broad_parent_significant %in% TRUE]
    }

    group_rows[[i]] <- data.frame(
      gene_set_group_id = sprintf("G%04d", i),
      exact_gene_set_shared_by_terms = nrow(x = group_df) > 1,
      term_count = nrow(x = group_df),
      gene_set_size = length(x = genes),
      geneID = paste(genes, collapse = "/"),
      representative_go_id = as.character(x = representative[[id_col]][[1]]),
      representative_description = as.character(x = representative[[term_col]][[1]]),
      representative_definition = representative_definition,
      best_pvalue = best_pvalue,
      best_adjusted_pvalue = best_padjust,
      ontologies = ontology_values,
      ontology_labels = ontology_label_values,
      go_ids = collapse_values(x = group_df[[id_col]]),
      go_terms = collapse_values(x = group_df[[term_col]]),
      go_definitions = go_definition_values,
      go_synonyms = go_synonym_values,
      immediate_parent_terms = immediate_parent_values,
      immediate_parent_summaries = immediate_parent_summary_values,
      significant_ancestor_terms = significant_ancestor_values,
      significant_descendant_terms = significant_descendant_values,
      bp_terms = collapse_values(x = bp_terms),
      mf_terms = collapse_values(x = mf_terms),
      cc_terms = collapse_values(x = cc_terms),
      most_specific_terms = collapse_values(x = most_specific_terms),
      broad_parent_terms = collapse_values(x = broad_parent_terms),
      stringsAsFactors = FALSE
    )
  }

  out <- do.call(what = rbind, args = group_rows)
  out <- out[
    order(
      -out$exact_gene_set_shared_by_terms,
      -out$term_count,
      out$best_adjusted_pvalue,
      out$best_pvalue,
      na.last = TRUE
    ),
    ,
    drop = FALSE
  ]
  rownames(x = out) <- NULL
  out
}

write_empty_gene_term_group_result <- function(file_path) {
  write_table_tsv(
    df = make_empty_gene_term_group_result(),
    file_path = file_path
  )
}

write_empty_result <- function(file_path, method) {
  if (identical(method, "clusterprofiler")) {
    empty_df <- data.frame(
      ID = character(),
      Description = character(),
      ontology = character(),
      ontology_label = character(),
      go_definition = character(),
      go_synonyms = character(),
      go_secondary_ids = character(),
      go_db_ontology = character(),
      go_immediate_parent_ids = character(),
      go_immediate_parent_terms = character(),
      go_immediate_parent_relations = character(),
      go_immediate_parent_summary = character(),
      go_ancestor_count = integer(),
      significant_descendant_count = integer(),
      has_significant_descendant = logical(),
      most_specific_significant = logical(),
      significant_ancestor_count = integer(),
      has_significant_ancestor = logical(),
      broad_parent_significant = logical(),
      significant_ancestor_terms = character(),
      significant_descendant_terms = character(),
      GeneRatio = character(),
      BgRatio = character(),
      pvalue = numeric(),
      p.adjust = numeric(),
      qvalue = numeric(),
      geneID = character(),
      Count = integer(),
      stringsAsFactors = FALSE
    )
  } else {
    empty_df <- data.frame(
      category = character(),
      ontology = character(),
      ontology_label = character(),
      go_definition = character(),
      go_synonyms = character(),
      go_secondary_ids = character(),
      go_db_ontology = character(),
      go_immediate_parent_ids = character(),
      go_immediate_parent_terms = character(),
      go_immediate_parent_relations = character(),
      go_immediate_parent_summary = character(),
      go_ancestor_count = integer(),
      significant_descendant_count = integer(),
      has_significant_descendant = logical(),
      most_specific_significant = logical(),
      significant_ancestor_count = integer(),
      has_significant_ancestor = logical(),
      broad_parent_significant = logical(),
      significant_ancestor_terms = character(),
      significant_descendant_terms = character(),
      over_represented_pvalue = numeric(),
      under_represented_pvalue = numeric(),
      numDEInCat = integer(),
      numInCat = integer(),
      over_represented_padj = numeric(),
      under_represented_padj = numeric(),
      go_term = character(),
      geneID = character(),
      Count = integer(),
      stringsAsFactors = FALSE
    )
  }

  utils::write.table(
    x = empty_df,
    file = file_path,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE,
    na = ""
  )
}

save_pdf_plot <- function(plot_fun, file_path) {
  grDevices::pdf(file = file_path, width = 10, height = 8)
  on.exit(expr = grDevices::dev.off(), add = TRUE)
  plot_fun()
}

make_clusterprofiler_plots <- function(result_df, plot_prefix, top_n_plot) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is not installed, so clusterProfiler PDF plots cannot be created.")
  }

  if (nrow(x = result_df) == 0) {
    return(invisible(NULL))
  }

  plot_df <- result_df
  plot_df <- plot_df[order(plot_df$p.adjust, plot_df$pvalue), , drop = FALSE]
  plot_df <- plot_df[!is.na(x = plot_df$p.adjust), , drop = FALSE]
  plot_df <- head(x = plot_df, n = top_n_plot)

  if (nrow(x = plot_df) == 0) {
    return(invisible(NULL))
  }

  plot_df$Description <- factor(
    x = plot_df$Description,
    levels = rev(x = unique(x = plot_df$Description))
  )
  plot_df$neg_log10_padj <- -log10(x = pmax(plot_df$p.adjust, .Machine$double.xmin))

  dot_file <- sprintf("%s_dotplot.pdf", plot_prefix)
  bar_file <- sprintf("%s_barplot.pdf", plot_prefix)

  save_pdf_plot(
    plot_fun = function() {
      p <- ggplot2::ggplot(
        data = plot_df,
        mapping = ggplot2::aes(x = neg_log10_padj, y = Description, size = Count)
      ) +
        ggplot2::geom_point() +
        ggplot2::labs(
          title = "clusterProfiler GO enrichment",
          x = expression(-log[10](adjusted~P)),
          y = "GO term",
          size = "Count"
        ) +
        ggplot2::theme_bw()
      print(p)
    },
    file_path = dot_file
  )

  save_pdf_plot(
    plot_fun = function() {
      p <- ggplot2::ggplot(
        data = plot_df,
        mapping = ggplot2::aes(x = neg_log10_padj, y = Description)
      ) +
        ggplot2::geom_col() +
        ggplot2::labs(
          title = "clusterProfiler GO enrichment",
          x = expression(-log[10](adjusted~P)),
          y = "GO term"
        ) +
        ggplot2::theme_bw()
      print(p)
    },
    file_path = bar_file
  )
}

make_goseq_plots <- function(result_df, pwf, plot_prefix, top_n_plot) {
  if (nrow(x = result_df) == 0) {
    return(invisible(NULL))
  }

  if (!is.null(x = pwf)) {
    pwf_file <- sprintf("%s_pwf.pdf", plot_prefix)
    save_pdf_plot(
      plot_fun = function() {
        goseq::plotPWF(pwf)
      },
      file_path = pwf_file
    )
  }

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is not installed, so goseq summary PDF plots cannot be created.")
  }

  plot_df <- result_df
  plot_df <- plot_df[!is.na(x = plot_df$over_represented_padj), , drop = FALSE]
  plot_df <- plot_df[order(plot_df$over_represented_padj, plot_df$over_represented_pvalue), , drop = FALSE]
  plot_df <- head(x = plot_df, n = top_n_plot)

  if (nrow(x = plot_df) == 0) {
    return(invisible(NULL))
  }

  plot_df$term_label <- ifelse(
    test = is.na(x = plot_df$go_term) | !nzchar(x = plot_df$go_term),
    yes = plot_df$category,
    no = plot_df$go_term
  )
  plot_df$term_label <- factor(
    x = plot_df$term_label,
    levels = rev(x = unique(x = plot_df$term_label))
  )
  plot_df$neg_log10_padj <- -log10(x = pmax(plot_df$over_represented_padj, .Machine$double.xmin))

  bar_file <- sprintf("%s_overrepresented_barplot.pdf", plot_prefix)

  save_pdf_plot(
    plot_fun = function() {
      p <- ggplot2::ggplot(
        data = plot_df,
        mapping = ggplot2::aes(x = neg_log10_padj, y = term_label)
      ) +
        ggplot2::geom_col() +
        ggplot2::labs(
          title = "goseq GO enrichment",
          x = expression(-log[10](adjusted~P)),
          y = "GO term"
        ) +
        ggplot2::theme_bw()
      print(p)
    },
    file_path = bar_file
  )
}

write_table_tsv <- function(df, file_path) {
  utils::write.table(
    x = df,
    file = file_path,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE,
    na = ""
  )
}


filter_significant_results <- function(result_df, method, fdr_cutoff) {
  if (nrow(x = result_df) == 0) {
    return(result_df)
  }

  if (identical(method, "clusterprofiler")) {
    if (!"p.adjust" %in% colnames(x = result_df)) {
      return(result_df[0, , drop = FALSE])
    }

    keep <- !is.na(x = result_df$p.adjust) & result_df$p.adjust <= fdr_cutoff
    return(result_df[keep, , drop = FALSE])
  }

  padj_cols <- intersect(
    x = c("over_represented_padj", "under_represented_padj"),
    y = colnames(x = result_df)
  )

  if (length(x = padj_cols) == 0) {
    return(result_df[0, , drop = FALSE])
  }

  keep_matrix <- do.call(
    what = cbind,
    args = lapply(
      X = padj_cols,
      FUN = function(col_name) {
        !is.na(x = result_df[[col_name]]) & result_df[[col_name]] <= fdr_cutoff
      }
    )
  )

  keep <- rowSums(x = keep_matrix) > 0
  result_df[keep, , drop = FALSE]
}

get_direction_dir_name <- function(direction_name) {
  direction_key <- tolower(x = trimws(x = direction_name))

  if (identical(direction_key, "up")) {
    return("upregulated")
  }

  if (identical(direction_key, "down")) {
    return("downregulated")
  }

  if (identical(direction_key, "all")) {
    return("all_significant")
  }

  safe_dir_name(x = direction_name)
}

create_output_paths <- function(out_dir, method_name, id_type, set_name, direction_name) {
  safe_set <- safe_dir_name(x = set_name)
  safe_direction <- safe_dir_name(x = direction_name)
  direction_dir_name <- get_direction_dir_name(direction_name = direction_name)

  method_dir <- file.path(
    out_dir,
    sprintf("%s_GO_enrichment_results", method_name)
  )
  direction_dir <- file.path(method_dir, direction_dir_name)
  significant_dir <- file.path(direction_dir, "significant")

  dir.create(path = direction_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(path = significant_dir, showWarnings = FALSE, recursive = TRUE)

  file_stub <- sprintf(
    "%s_%s_%s_%s",
    method_name,
    id_type,
    safe_set,
    safe_direction
  )

  list(
    method_dir = method_dir,
    direction_dir = direction_dir,
    significant_dir = significant_dir,
    full_result_file = file.path(direction_dir, sprintf("%s_go_enrichment.tsv", file_stub)),
    full_plot_prefix = file.path(direction_dir, file_stub),
    significant_result_file = file.path(
      significant_dir,
      sprintf("%s_significant_go_enrichment.tsv", file_stub)
    ),
    specific_result_file = file.path(
      significant_dir,
      sprintf("%s_significant_most_specific_go_enrichment.tsv", file_stub)
    ),
    broad_result_file = file.path(
      significant_dir,
      sprintf("%s_significant_broad_parent_go_enrichment.tsv", file_stub)
    ),
    gene_group_result_file = file.path(
      significant_dir,
      sprintf("%s_significant_gene_set_term_groups.tsv", file_stub)
    ),
    significant_plot_prefix = file.path(significant_dir, sprintf("%s_significant", file_stub)),
    specific_plot_prefix = file.path(significant_dir, sprintf("%s_significant_most_specific", file_stub)),
    broad_plot_prefix = file.path(significant_dir, sprintf("%s_significant_broad_parent", file_stub))
  )
}

check_packages(methods = args$methods)
dir.create(path = args$out_dir, showWarnings = FALSE, recursive = TRUE)

warning_messages <- character()
info_messages <- character()

log_info <- function(msg) {
  timestamp <- format(x = Sys.time(), format = "%Y-%m-%d %H:%M:%S")
  line <- sprintf("[%s] INFO: %s", timestamp, msg)
  info_messages <<- c(info_messages, line)
  message(line)
}

log_warning <- function(msg) {
  timestamp <- format(x = Sys.time(), format = "%Y-%m-%d %H:%M:%S")
  line <- sprintf("[%s] WARNING: %s", timestamp, msg)
  warning_messages <<- c(warning_messages, line)
  message(line)
}

log_info("Reading annotation file.")
annotation_info <- parse_annotation(
  annotation_file = args$annotation_file,
  id_type = args$id_type,
  ontology = args$ontology
)

if (!is.null(x = args$de_file)) {
  log_info("Reading differential expression file.")
  de_df <- read_tabular_file(file_path = args$de_file)
} else {
  de_df <- NULL
}

background_ids <- get_background_ids(
  annotation_info = annotation_info,
  de_df = de_df,
  background_source = args$background_source,
  background_file = args$background_file,
  expression_file = args$expression_file,
  expression_min_tpm = args$expression_min_tpm,
  expression_min_samples = args$expression_min_samples
)


annot_universe <- unique(
  trimws(x = as.character(x = annotation_info$annotation_df[[annotation_info$id_col]]))
)

log_info(sprintf("Background source: %s", args$background_source))
if (!is.null(x = args$expression_file)) {
  log_info(sprintf("Expression file: %s", args$expression_file))
}
log_info(sprintf("Background IDs before intersect: %s", length(x = background_ids)))
log_info(sprintf(
  "First 5 background IDs: %s",
  paste(utils::head(x = background_ids, n = 5), collapse = ", ")
))
log_info(sprintf("Annotation ID column: %s", annotation_info$id_col))
log_info(sprintf("Annotation IDs before intersect: %s", length(x = annot_universe)))
log_info(sprintf(
  "First 5 annotation IDs: %s",
  paste(utils::head(x = annot_universe, n = 5), collapse = ", ")
))
log_info(sprintf(
  "Overlap before intersect: %s",
  length(x = intersect(
    x = trimws(x = as.character(x = background_ids)),
    y = annot_universe
  ))
))

background_ids <- intersect(
  x = trimws(x = as.character(x = background_ids)),
  y = annot_universe
)

if (length(x = background_ids) == 0) {
  stop(
    paste(
      "No background IDs remained after intersecting with the annotation file.",
      "Please inspect the debug lines above for ID-format mismatches."
    ),
    call. = FALSE
  )
}

log_info(sprintf("Background universe size: %s", length(x = background_ids)))

term2gene <- annotation_info$term2gene[
  annotation_info$term2gene$feature_id %in% background_ids,
  ,
  drop = FALSE
]
term2name <- annotation_info$term2name[
  annotation_info$term2name$go_id %in% term2gene$go_id,
  ,
  drop = FALSE
]

if (nrow(x = term2gene) == 0) {
  stop("No GO mappings remained after filtering to the background universe.", call. = FALSE)
}

input_sets <- list()

if (!is.null(x = args$de_file)) {
  input_sets <- prepare_de_sets(
    de_df = de_df,
    padj_cutoff = args$padj_cutoff,
    abs_log2fc_cutoff = args$abs_log2fc_cutoff,
    split_direction = args$split_direction
  )
} else {
  log_info("Preparing input ID list sets.")
  input_sets <- prepare_id_list_sets(id_list_file = args$id_list_file)
}

length_vector <- NULL
if ("goseq" %in% args$methods) {
  length_vector <- get_length_vector(
    background_ids = background_ids,
    fasta_file = args$fasta_file,
    annotation_info = annotation_info,
    id_type = args$id_type
  )
}

summary_rows <- list()

for (method_name in args$methods) {
  dir.create(
    path = file.path(args$out_dir, sprintf("%s_GO_enrichment_results", method_name)),
    showWarnings = FALSE,
    recursive = TRUE
  )
}

for (set_name in names(x = input_sets)) {
  log_info(sprintf("Processing set: %s", set_name))

  for (direction_name in names(x = input_sets[[set_name]])) {
    raw_ids <- input_sets[[set_name]][[direction_name]]
    sig_ids <- intersect(x = raw_ids, y = background_ids)

    if (length(x = raw_ids) == 0) {
      log_warning(sprintf(
        "Input set '%s' and direction '%s' contained no IDs.",
        set_name,
        direction_name
      ))
    }

    if (length(x = sig_ids) < length(x = raw_ids)) {
      log_warning(sprintf(
        paste(
          "For set '%s' and direction '%s', %s of %s input IDs matched the",
          "selected background and annotation universe."
        ),
        set_name,
        direction_name,
        length(x = sig_ids),
        length(x = raw_ids)
      ))
    }

    for (method_name in args$methods) {
      output_paths <- create_output_paths(
        out_dir = args$out_dir,
        method_name = method_name,
        id_type = args$id_type,
        set_name = set_name,
        direction_name = direction_name
      )

      full_result_file <- output_paths$full_result_file
      significant_result_file <- output_paths$significant_result_file
      specific_result_file <- output_paths$specific_result_file
      broad_result_file <- output_paths$broad_result_file
      gene_group_result_file <- output_paths$gene_group_result_file
      full_plot_prefix <- output_paths$full_plot_prefix
      significant_plot_prefix <- output_paths$significant_plot_prefix
      specific_plot_prefix <- output_paths$specific_plot_prefix
      broad_plot_prefix <- output_paths$broad_plot_prefix

      status <- "success"
      note <- ""
      result_count <- 0
      significant_result_count <- 0
      gene_term_group_count <- 0

      if (length(x = sig_ids) == 0) {
        write_empty_result(file_path = full_result_file, method = method_name)
        write_empty_result(file_path = significant_result_file, method = method_name)
        write_empty_result(file_path = specific_result_file, method = method_name)
        write_empty_result(file_path = broad_result_file, method = method_name)
        write_empty_gene_term_group_result(file_path = gene_group_result_file)
        status <- "empty"
        note <- "No input IDs remained after filtering to the background universe."

        summary_rows[[length(summary_rows) + 1]] <- data.frame(
          set_name = set_name,
          direction = direction_name,
          method = method_name,
          input_id_count = length(x = raw_ids),
          analysed_id_count = 0,
          background_count = length(x = background_ids),
          result_count = 0,
          significant_result_count = 0,
          gene_term_group_count = 0,
          fdr_cutoff = args$fdr_cutoff,
          status = status,
          note = note,
          full_result_file = full_result_file,
          significant_result_file = significant_result_file,
          specific_result_file = specific_result_file,
          broad_result_file = broad_result_file,
          gene_group_result_file = gene_group_result_file,
          direction_dir = output_paths$direction_dir,
          significant_dir = output_paths$significant_dir,
          stringsAsFactors = FALSE
        )
        next
      }

      analysis_output <- tryCatch(
        expr = {
          if (identical(method_name, "clusterprofiler")) {
            run_clusterprofiler(
              sig_ids = sig_ids,
              background_ids = background_ids,
              term2gene = term2gene,
              term2name = term2name,
              min_gs_size = args$min_gs_size,
              max_gs_size = args$max_gs_size
            )
          } else {
            run_goseq(
              sig_ids = sig_ids,
              background_ids = background_ids,
              term2gene = term2gene,
              term2name = term2name,
              length_vector = length_vector
            )
          }
        },
        error = function(e) {
          log_warning(sprintf(
            "Analysis failed for set '%s', direction '%s', method '%s': %s",
            set_name,
            direction_name,
            method_name,
            conditionMessage(e)
          ))
          NULL
        }
      )

      if (is.null(x = analysis_output)) {
        write_empty_result(file_path = full_result_file, method = method_name)
        write_empty_result(file_path = significant_result_file, method = method_name)
        write_empty_result(file_path = specific_result_file, method = method_name)
        write_empty_result(file_path = broad_result_file, method = method_name)
        write_empty_gene_term_group_result(file_path = gene_group_result_file)
        status <- "analysis_failed"
        note <- "Analysis failed. See warnings log."
      } else {
        result_df <- analysis_output$result_df
        result_df <- add_go_context_columns(
          result_df = result_df,
          method = method_name,
          term2name = term2name
        )
        result_count <- nrow(x = result_df)
        significant_df <- filter_significant_results(
          result_df = result_df,
          method = method_name,
          fdr_cutoff = args$fdr_cutoff
        )
        significant_df <- mark_most_specific_significant(
          significant_df = significant_df,
          method = method_name
        )
        specific_df <- filter_most_specific_significant(
          significant_df = significant_df,
          method = method_name
        )
        broad_df <- filter_broad_parent_significant(
          significant_df = significant_df,
          method = method_name
        )
        gene_term_group_df <- make_gene_set_term_groups(
          significant_df = significant_df,
          method = method_name
        )
        significant_result_count <- nrow(x = significant_df)
        gene_term_group_count <- nrow(x = gene_term_group_df)

        if (nrow(x = result_df) == 0) {
          write_empty_result(file_path = full_result_file, method = method_name)
          status <- "no_terms"
          note <- "No enriched GO terms were returned."
        } else {
          write_table_tsv(df = result_df, file_path = full_result_file)
        }

        if (nrow(x = significant_df) == 0) {
          write_empty_result(file_path = significant_result_file, method = method_name)

          if (identical(status, "success")) {
            status <- "no_significant_terms"
            note <- sprintf(
              "No GO terms passed the FDR cutoff of %s.",
              args$fdr_cutoff
            )
          }
        } else {
          write_table_tsv(df = significant_df, file_path = significant_result_file)
        }

        if (args$write_specific_terms) {
          if (exists(x = "specific_df") && nrow(x = specific_df) > 0) {
            write_table_tsv(df = specific_df, file_path = specific_result_file)
          } else {
            write_empty_result(file_path = specific_result_file, method = method_name)
          }
        }

        if (args$write_broad_terms) {
          if (exists(x = "broad_df") && nrow(x = broad_df) > 0) {
            write_table_tsv(df = broad_df, file_path = broad_result_file)
          } else {
            write_empty_result(file_path = broad_result_file, method = method_name)
          }
        }

        if (args$write_gene_term_groups) {
          if (exists(x = "gene_term_group_df") && nrow(x = gene_term_group_df) > 0) {
            write_table_tsv(df = gene_term_group_df, file_path = gene_group_result_file)
          } else {
            write_empty_gene_term_group_result(file_path = gene_group_result_file)
          }
        }

        if (args$make_plots) {
          tryCatch(
            expr = {
              if (identical(method_name, "clusterprofiler")) {
                make_clusterprofiler_plots(
                  result_df = result_df,
                  plot_prefix = full_plot_prefix,
                  top_n_plot = args$top_n_plot
                )
              } else {
                make_goseq_plots(
                  result_df = result_df,
                  pwf = analysis_output$pwf,
                  plot_prefix = full_plot_prefix,
                  top_n_plot = args$top_n_plot
                )
              }
            },
            error = function(e) {
              log_warning(sprintf(
                "Plotting failed for set '%s', direction '%s', method '%s': %s",
                set_name,
                direction_name,
                method_name,
                conditionMessage(e)
              ))
              invisible(NULL)
            }
          )

          if (nrow(x = significant_df) > 0) {
            tryCatch(
              expr = {
                if (identical(method_name, "clusterprofiler")) {
                  make_clusterprofiler_plots(
                    result_df = significant_df,
                    plot_prefix = significant_plot_prefix,
                    top_n_plot = args$top_n_plot
                  )
                } else {
                  make_goseq_plots(
                    result_df = significant_df,
                    pwf = analysis_output$pwf,
                    plot_prefix = significant_plot_prefix,
                    top_n_plot = args$top_n_plot
                  )
                }
              },
              error = function(e) {
                log_warning(sprintf(
                  paste(
                    "Significant-result plotting failed for set '%s', direction '%s',",
                    "method '%s': %s"
                  ),
                  set_name,
                  direction_name,
                  method_name,
                  conditionMessage(e)
                ))
                invisible(NULL)
              }
            )
          }
        }
      }

      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        set_name = set_name,
        direction = direction_name,
        method = method_name,
        input_id_count = length(x = raw_ids),
        analysed_id_count = length(x = sig_ids),
        background_count = length(x = background_ids),
        result_count = result_count,
        significant_result_count = significant_result_count,
        gene_term_group_count = gene_term_group_count,
        fdr_cutoff = args$fdr_cutoff,
        status = status,
        note = note,
        full_result_file = full_result_file,
        significant_result_file = significant_result_file,
        specific_result_file = specific_result_file,
        broad_result_file = broad_result_file,
        gene_group_result_file = gene_group_result_file,
        direction_dir = output_paths$direction_dir,
        significant_dir = output_paths$significant_dir,
        stringsAsFactors = FALSE
      )
    }
  }
}

summary_df <- do.call(what = rbind, args = summary_rows)
summary_file <- file.path(path = args$out_dir, "go_enrichment_run_summary.tsv")
write_table_tsv(df = summary_df, file_path = summary_file)

args_file <- file.path(path = args$out_dir, "run_arguments.tsv")
args_df <- data.frame(
  argument = names(x = args),
  value = vapply(
    X = args,
    FUN = function(x) {
      if (length(x = x) == 0 || is.null(x = x)) {
        return("")
      }
      paste(x, collapse = ",")
    },
    FUN.VALUE = character(length = 1)
  ),
  stringsAsFactors = FALSE
)
write_table_tsv(df = args_df, file_path = args_file)

session_file <- file.path(path = args$out_dir, "session_info.txt")
writeLines(text = capture.output(utils::sessionInfo()), con = session_file)

info_file <- file.path(path = args$out_dir, "run_log.txt")
warning_file <- file.path(path = args$out_dir, "warnings_log.txt")

log_info("Finished.")
log_info(sprintf("Summary written to: %s", summary_file))

if (length(x = info_messages) == 0) {
  writeLines(text = "No info messages recorded.", con = info_file)
} else {
  writeLines(text = info_messages, con = info_file)
}

if (length(x = warning_messages) == 0) {
  writeLines(text = "No warnings recorded.", con = warning_file)
} else {
  writeLines(text = warning_messages, con = warning_file)
}
