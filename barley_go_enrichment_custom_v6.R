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
  seq_lengths <- width(x = seqs)
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
    "go_ancestor_count",
    "significant_descendant_count",
    "has_significant_descendant",
    "most_specific_significant",
    "significant_ancestor_count",
    "has_significant_ancestor",
    "broad_parent_significant"
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

  result_df
}

mark_go_hierarchy_significance <- function(significant_df, method) {
  if (nrow(x = significant_df) == 0) {
    return(significant_df)
  }

  id_col <- if (identical(method, "clusterprofiler")) "ID" else "category"

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
    descendant_count[significant_ancestors] <- descendant_count[significant_ancestors] + 1L
  }

  significant_df$significant_descendant_count <- as.integer(x = descendant_count[go_ids])
  significant_df$has_significant_descendant <- significant_df$significant_descendant_count > 0L
  significant_df$most_specific_significant <- !significant_df$has_significant_descendant
  significant_df$significant_ancestor_count <- as.integer(x = ancestor_count[go_ids])
  significant_df$has_significant_ancestor <- significant_df$significant_ancestor_count > 0L
  significant_df$broad_parent_significant <- !significant_df$has_significant_ancestor

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

write_empty_result <- function(file_path, method) {
  if (identical(method, "clusterprofiler")) {
    empty_df <- data.frame(
      ID = character(),
      Description = character(),
      ontology = character(),
      ontology_label = character(),
      go_ancestor_count = integer(),
      significant_descendant_count = integer(),
      has_significant_descendant = logical(),
      most_specific_significant = logical(),
      significant_ancestor_count = integer(),
      has_significant_ancestor = logical(),
      broad_parent_significant = logical(),
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
      go_ancestor_count = integer(),
      significant_descendant_count = integer(),
      has_significant_descendant = logical(),
      most_specific_significant = logical(),
      significant_ancestor_count = integer(),
      has_significant_ancestor = logical(),
      broad_parent_significant = logical(),
      over_represented_pvalue = numeric(),
      under_represented_pvalue = numeric(),
      numDEInCat = integer(),
      numInCat = integer(),
      over_represented_padj = numeric(),
      under_represented_padj = numeric(),
      go_term = character(),
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
      full_plot_prefix <- output_paths$full_plot_prefix
      significant_plot_prefix <- output_paths$significant_plot_prefix
      specific_plot_prefix <- output_paths$specific_plot_prefix
      broad_plot_prefix <- output_paths$broad_plot_prefix

      status <- "success"
      note <- ""
      result_count <- 0
      significant_result_count <- 0

      if (length(x = sig_ids) == 0) {
        write_empty_result(file_path = full_result_file, method = method_name)
        write_empty_result(file_path = significant_result_file, method = method_name)
        write_empty_result(file_path = specific_result_file, method = method_name)
        write_empty_result(file_path = broad_result_file, method = method_name)
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
          fdr_cutoff = args$fdr_cutoff,
          status = status,
          note = note,
          full_result_file = full_result_file,
          significant_result_file = significant_result_file,
          specific_result_file = specific_result_file,
          broad_result_file = broad_result_file,
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
        significant_result_count <- nrow(x = significant_df)

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
        fdr_cutoff = args$fdr_cutoff,
        status = status,
        note = note,
        full_result_file = full_result_file,
        significant_result_file = significant_result_file,
        specific_result_file = specific_result_file,
        broad_result_file = broad_result_file,
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
