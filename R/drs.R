#' Prepare an expression matrix for DRS scoring
#'
#' `data_pre()` reads a gene-by-sample expression matrix from a file or
#' standardizes an in-memory matrix/data frame. By default, rows are genes
#' and columns are samples. Duplicate gene symbols are merged by row-wise mean.
#'
#' @param x A file path to a csv/tsv/txt matrix, or an in-memory matrix/data
#'   frame.
#' @param gene_col For file or data frame inputs where gene symbols are stored
#'   in a column, the 1-based index of that column. Defaults to `1`.
#' @param sep Optional delimiter override. By default the delimiter is inferred
#'   from the file extension or file contents.
#'
#' @return A numeric matrix with gene symbols as row names and samples as
#'   columns.
#' @export
data_pre <- function(x, gene_col = 1, sep = NULL) {
  .prepare_expression_input(x = x, gene_col = gene_col, sep = sep)
}

#' Load and validate the bundled DRS coefficient model
#'
#' `load_drs_model()` reads `drs_coefficients.csv` from `inst/extdata/` and
#' validates that the coefficient table contains the required `Gene` and
#' `Coefficient` columns for the complete 20-gene DRS model.
#'
#' @param file Optional path to a coefficient csv file. When `NULL`, the
#'   bundled `drs_coefficients.csv` file is used.
#'
#' @return A data frame with columns `Gene` and `Coefficient`.
#' @export
load_drs_model <- function(file = NULL) {
  if (is.null(file)) {
    file <- .locate_extdata_file("drs_coefficients.csv")
  }

  if (!file.exists(file)) {
    stop("Coefficient file does not exist: ", file, call. = FALSE)
  }

  model_df <- utils::read.csv(
    file = file,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  .validate_model_dataframe(model_df)
}

#' Predict weighted DRS scores for each sample
#'
#' `predict_drs()` applies the bundled 20-gene coefficient model to an
#' expression matrix. Each model gene is z-scored across the input cohort
#' before the weighted sum is calculated.
#'
#' @param input_data A prepared numeric matrix, a data frame, or a path to an
#'   expression matrix file.
#' @param model Optional coefficient data frame returned by `load_drs_model()`.
#'   When `NULL`, the bundled model is used.
#' @param group_cutoff Numeric cutoff used to define `DRS_group`. Scores greater
#'   than or equal to the cutoff are labeled `"High"`; lower scores are labeled
#'   `"Low"`. Defaults to `0`.
#'
#' @return A data frame with columns `Sample`, `DRS_score`, `DRS_group`,
#'   `n_matched_genes`, `n_missing_genes`, and `missing_genes`.
#' @export
predict_drs <- function(input_data, model = NULL, group_cutoff = 0) {
  if (is.null(model)) {
    model <- load_drs_model()
  } else {
    model <- .validate_model_dataframe(model)
  }

  expression_mat <- .prepare_expression_input(input_data)
  .validate_group_cutoff(group_cutoff)

  if (ncol(expression_mat) < 2L) {
    stop("At least two samples are required to compute gene-wise z-scores.", call. = FALSE)
  }

  matched_genes <- intersect(model$Gene, rownames(expression_mat))
  missing_genes <- setdiff(model$Gene, rownames(expression_mat))

  if (length(missing_genes) > 0L) {
    stop(
      "Missing required model genes: ",
      paste(missing_genes, collapse = ", "),
      call. = FALSE
    )
  }

  model_expression <- expression_mat[model$Gene, , drop = FALSE]

  if (anyNA(model_expression)) {
    na_genes <- rownames(model_expression)[apply(is.na(model_expression), 1L, any)]
    stop(
      "Missing expression values were detected for model genes: ",
      paste(na_genes, collapse = ", "),
      call. = FALSE
    )
  }

  gene_sd <- apply(model_expression, 1L, stats::sd)
  zero_variance_genes <- names(gene_sd)[is.na(gene_sd) | gene_sd == 0]

  if (length(zero_variance_genes) > 0L) {
    stop(
      "Zero-variance model genes detected: ",
      paste(zero_variance_genes, collapse = ", "),
      call. = FALSE
    )
  }

  gene_mean <- rowMeans(model_expression)
  z_mat <- sweep(model_expression, 1L, gene_mean, FUN = "-")
  z_mat <- sweep(z_mat, 1L, gene_sd, FUN = "/")
  weighted_mat <- sweep(z_mat, 1L, model$Coefficient, FUN = "*")
  drs_scores <- colSums(weighted_mat)

  sample_names <- colnames(model_expression)
  if (is.null(sample_names) || any(!nzchar(sample_names))) {
    sample_names <- paste0("Sample", seq_len(ncol(model_expression)))
  }

  data.frame(
    Sample = sample_names,
    DRS_score = as.numeric(drs_scores),
    DRS_group = ifelse(drs_scores >= group_cutoff, "High", "Low"),
    n_matched_genes = length(matched_genes),
    n_missing_genes = length(missing_genes),
    missing_genes = if (length(missing_genes) > 0L) {
      paste(missing_genes, collapse = ",")
    } else {
      ""
    },
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

#' Convenience wrapper for DRS prediction
#'
#' `drs_pre()` is a lightweight wrapper around [predict_drs()] that keeps the
#' user workflow concise in examples and scripts.
#'
#' @param input_data A prepared numeric matrix, a data frame, or a path to an
#'   expression matrix file.
#' @param model Optional coefficient data frame returned by `load_drs_model()`.
#'   When `NULL`, the bundled model is used.
#' @param group_cutoff Numeric cutoff used to define `DRS_group`. Defaults to
#'   `0`.
#'
#' @return A data frame with DRS prediction results.
#' @export
drs_pre <- function(input_data, model = NULL, group_cutoff = 0) {
  predict_drs(input_data = input_data, model = model, group_cutoff = group_cutoff)
}

#' Plot the DRS score distribution
#'
#' `plot_drs_distribution()` draws a simple base-R bar plot of sample-level DRS
#' scores ordered from low to high. The function can accept a prediction result
#' from [predict_drs()] or any input supported by [drs_pre()].
#'
#' @param x A prediction data frame returned by [predict_drs()], or any input
#'   accepted by [drs_pre()].
#' @param model Optional coefficient data frame returned by `load_drs_model()`
#'   when `x` is not already a prediction result.
#' @param group_cutoff Numeric cutoff used to define `DRS_group` when `x` is not
#'   already a prediction result. Defaults to `0`.
#' @param main Plot title.
#' @param ylab Y-axis label.
#' @param las Axis label orientation passed to [graphics::barplot()].
#' @param ... Additional arguments passed to [graphics::barplot()].
#'
#' @return Invisibly returns the plotted prediction data frame.
#' @export
plot_drs_distribution <- function(
  x,
  model = NULL,
  group_cutoff = 0,
  main = "DRS score distribution",
  ylab = "DRS score",
  las = 2,
  ...
) {
  if (is.data.frame(x) && all(c("Sample", "DRS_score", "DRS_group") %in% names(x))) {
    result <- x
  } else {
    result <- drs_pre(input_data = x, model = model, group_cutoff = group_cutoff)
  }

  ordered_result <- result[order(result$DRS_score), , drop = FALSE]
  bar_col <- ifelse(ordered_result$DRS_group == "High", "#D55E00", "#0072B2")

  graphics::barplot(
    height = ordered_result$DRS_score,
    names.arg = ordered_result$Sample,
    col = bar_col,
    border = NA,
    las = las,
    ylab = ylab,
    main = main,
    ...
  )
  graphics::abline(h = group_cutoff, lty = 2, lwd = 2, col = "gray30")
  graphics::legend(
    "topright",
    legend = c("High", "Low"),
    fill = c("#D55E00", "#0072B2"),
    border = NA,
    bty = "n"
  )

  invisible(ordered_result)
}

.drs_model_template <- function() {
  data.frame(
    Gene = c(
      "RAC1", "CYFIP1", "WASF2", "BRK1", "ABI1",
      "ACTR2", "ARPC2", "FLNA", "FLNB", "ACTN4",
      "IQGAP1", "MYH10", "OSTC", "ACSL4", "DDO",
      "PDLIM1", "DSTN", "VCL", "CD44", "YTHDC1"
    ),
    Coefficient = c(
      0.0477714594953878, 0.0845378169964754, 0.228860923807596,
      0.111178465923049, -0.690017908034411, 0.0608166918319247,
      0.030483537785947, 0.159431512080148, -0.280178147690666,
      0.3770931565415, 0.186772792281645, 0.0234811083291907,
      0.404433914021094, 0.345430775834607, 0.274010134974565,
      0.0206130110104642, -0.232575375529517, -0.25377594193562,
      0.0305987238207958, -0.376141806552824
    ),
    stringsAsFactors = FALSE
  )
}

.validate_model_dataframe <- function(model_df) {
  if (!is.data.frame(model_df)) {
    stop("The model must be supplied as a data frame.", call. = FALSE)
  }

  required_cols <- c("Gene", "Coefficient")
  missing_cols <- setdiff(required_cols, colnames(model_df))
  if (length(missing_cols) > 0L) {
    stop(
      "The coefficient table must contain columns: Gene and Coefficient. Missing: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  model_df <- model_df[, required_cols, drop = FALSE]
  model_df$Gene <- toupper(trimws(as.character(model_df$Gene)))

  suppressWarnings(coef_numeric <- as.numeric(model_df$Coefficient))
  invalid_coef <- is.na(coef_numeric)
  if (any(invalid_coef)) {
    stop("All coefficient values must be numeric.", call. = FALSE)
  }
  model_df$Coefficient <- coef_numeric

  if (anyDuplicated(model_df$Gene) > 0L) {
    duplicated_genes <- unique(model_df$Gene[duplicated(model_df$Gene)])
    stop(
      "Duplicate genes were found in the coefficient table: ",
      paste(duplicated_genes, collapse = ", "),
      call. = FALSE
    )
  }

  template <- .drs_model_template()
  missing_genes <- setdiff(template$Gene, model_df$Gene)
  extra_genes <- setdiff(model_df$Gene, template$Gene)

  if (length(missing_genes) > 0L || length(extra_genes) > 0L) {
    problems <- character()
    if (length(missing_genes) > 0L) {
      problems <- c(problems, paste0("missing genes: ", paste(missing_genes, collapse = ", ")))
    }
    if (length(extra_genes) > 0L) {
      problems <- c(problems, paste0("unexpected genes: ", paste(extra_genes, collapse = ", ")))
    }
    stop(
      "The coefficient table does not match the required 20-gene DRS model (",
      paste(problems, collapse = "; "),
      ").",
      call. = FALSE
    )
  }

  ordered_idx <- match(template$Gene, model_df$Gene)
  model_df[ordered_idx, , drop = FALSE]
}

.prepare_expression_input <- function(x, gene_col = 1, sep = NULL) {
  if (is.character(x) && length(x) == 1L) {
    if (!file.exists(x)) {
      stop("Input file does not exist: ", x, call. = FALSE)
    }
    return(.read_expression_file(path = x, gene_col = gene_col, sep = sep))
  }

  if (is.matrix(x)) {
    return(.standardize_matrix_input(x))
  }

  if (is.data.frame(x)) {
    return(.standardize_data_frame_input(x, gene_col = gene_col))
  }

  stop(
    "Input must be a file path, matrix, or data frame.",
    call. = FALSE
  )
}

.read_expression_file <- function(path, gene_col = 1, sep = NULL) {
  if (!file.exists(path)) {
    stop("Input file does not exist: ", path, call. = FALSE)
  }

  if (is.null(sep)) {
    sep <- .infer_separator(path)
  }

  expression_df <- utils::read.table(
    file = path,
    header = TRUE,
    sep = sep,
    check.names = FALSE,
    quote = "\"",
    comment.char = "",
    stringsAsFactors = FALSE
  )

  .standardize_data_frame_input(expression_df, gene_col = gene_col)
}

.standardize_data_frame_input <- function(x, gene_col = 1) {
  if (!is.numeric(gene_col) || length(gene_col) != 1L || gene_col < 1L || gene_col > ncol(x)) {
    stop("`gene_col` must identify a valid gene column.", call. = FALSE)
  }

  gene_symbols <- x[[gene_col]]
  expression_df <- x[, -gene_col, drop = FALSE]

  if (ncol(expression_df) < 2L) {
    stop("At least two sample columns are required.", call. = FALSE)
  }

  .matrix_from_gene_column(expression_df = expression_df, gene_symbols = gene_symbols)
}

.standardize_matrix_input <- function(x) {
  if (is.null(rownames(x))) {
    stop("Matrix input must contain gene symbols as row names.", call. = FALSE)
  }

  if (ncol(x) < 2L) {
    stop("At least two sample columns are required.", call. = FALSE)
  }

  raw_vec <- as.vector(x)
  suppressWarnings(num_vec <- as.numeric(raw_vec))
  invalid_mask <- is.na(num_vec) &
    !is.na(raw_vec) &
    !(trimws(as.character(raw_vec)) %in% c("", "NA", "NaN"))

  if (any(invalid_mask)) {
    bad_value <- unique(as.character(raw_vec[invalid_mask]))[1L]
    stop("Non-numeric expression values detected: ", bad_value, call. = FALSE)
  }

  numeric_mat <- matrix(
    num_vec,
    nrow = nrow(x),
    ncol = ncol(x),
    dimnames = dimnames(x)
  )

  .postprocess_expression_matrix(numeric_mat)
}

.matrix_from_gene_column <- function(expression_df, gene_symbols) {
  raw_mat <- as.matrix(expression_df)
  raw_vec <- as.vector(raw_mat)
  suppressWarnings(num_vec <- as.numeric(raw_vec))
  invalid_mask <- is.na(num_vec) &
    !is.na(raw_vec) &
    !(trimws(as.character(raw_vec)) %in% c("", "NA", "NaN"))

  if (any(invalid_mask)) {
    bad_value <- unique(as.character(raw_vec[invalid_mask]))[1L]
    stop("Non-numeric expression values detected: ", bad_value, call. = FALSE)
  }

  numeric_mat <- matrix(
    num_vec,
    nrow = nrow(raw_mat),
    ncol = ncol(raw_mat),
    dimnames = list(gene_symbols, colnames(expression_df))
  )

  .postprocess_expression_matrix(numeric_mat)
}

.postprocess_expression_matrix <- function(x) {
  rownames(x) <- toupper(trimws(as.character(rownames(x))))

  if (any(!nzchar(rownames(x)))) {
    stop("Gene symbols must not be empty.", call. = FALSE)
  }

  if (is.null(colnames(x)) || any(!nzchar(colnames(x)))) {
    colnames(x) <- paste0("Sample", seq_len(ncol(x)))
  }

  if (anyDuplicated(colnames(x)) > 0L) {
    duplicated_samples <- unique(colnames(x)[duplicated(colnames(x))])
    stop(
      "Duplicate sample names were detected: ",
      paste(duplicated_samples, collapse = ", "),
      call. = FALSE
    )
  }

  .collapse_duplicate_genes(x)
}

.collapse_duplicate_genes <- function(x) {
  if (anyDuplicated(rownames(x)) == 0L) {
    return(x)
  }

  unique_genes <- unique(rownames(x))
  collapsed <- t(vapply(
    unique_genes,
    FUN = function(gene) {
      cols <- x[rownames(x) == gene, , drop = FALSE]
      colMeans(cols, na.rm = TRUE)
    },
    FUN.VALUE = numeric(ncol(x))
  ))

  rownames(collapsed) <- unique_genes
  colnames(collapsed) <- colnames(x)
  collapsed
}

.infer_separator <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (identical(ext, "csv")) {
    return(",")
  }
  if (ext %in% c("tsv", "tab")) {
    return("\t")
  }

  lines <- readLines(path, n = 5L, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]

  if (length(lines) == 0L) {
    stop("Input file is empty: ", path, call. = FALSE)
  }

  probe <- lines[1L]
  delimiter_counts <- c(
    comma = lengths(regmatches(probe, gregexpr(",", probe, fixed = TRUE))),
    tab = lengths(regmatches(probe, gregexpr("\t", probe, fixed = TRUE))),
    semicolon = lengths(regmatches(probe, gregexpr(";", probe, fixed = TRUE)))
  )

  delimiter_counts[delimiter_counts < 0L] <- 0L

  if (max(delimiter_counts) == 0L) {
    ""
  } else if (which.max(delimiter_counts) == 1L) {
    ","
  } else if (which.max(delimiter_counts) == 2L) {
    "\t"
  } else {
    ";"
  }
}

.validate_group_cutoff <- function(group_cutoff) {
  if (!is.numeric(group_cutoff) || length(group_cutoff) != 1L || is.na(group_cutoff)) {
    stop("`group_cutoff` must be a single numeric value.", call. = FALSE)
  }
}

.locate_extdata_file <- function(file_name) {
  installed_path <- system.file("extdata", file_name, package = "DisulfRisk")
  if (nzchar(installed_path) && file.exists(installed_path)) {
    return(installed_path)
  }

  root <- .find_package_root()
  if (!is.null(root)) {
    source_path <- file.path(root, "inst", "extdata", file_name)
    if (file.exists(source_path)) {
      return(source_path)
    }
  }

  stop("Could not locate extdata file: ", file_name, call. = FALSE)
}

.find_package_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = FALSE)

  repeat {
    if (file.exists(file.path(current, "DESCRIPTION"))) {
      return(current)
    }

    parent <- dirname(current)
    if (identical(parent, current)) {
      return(NULL)
    }
    current <- parent
  }
}
