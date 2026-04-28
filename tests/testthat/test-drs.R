example_expression_path <- function() {
  path <- system.file("extdata", "example_expression.csv", package = "DisulfRisk")
  if (!nzchar(path)) {
    path <- testthat::test_path("../../inst/extdata/example_expression.csv")
  }
  path
}

test_that("load_drs_model validates the bundled 20-gene model", {
  model <- load_drs_model()

  expect_s3_class(model, "data.frame")
  expect_identical(colnames(model), c("Gene", "Coefficient"))
  expect_equal(nrow(model), 20L)
  expect_true(all(c("RAC1", "YTHDC1") %in% model$Gene))
})

test_that("data_pre reads the example expression matrix", {
  path <- example_expression_path()
  input_data <- data_pre(path)

  expect_true(is.matrix(input_data))
  expect_equal(dim(input_data), c(24L, 6L))
  expect_true(all(c("RAC1", "ACTB") %in% rownames(input_data)))
  expect_true(all(colnames(input_data) == paste0("Sample_0", 1:6)))
})

test_that("data_pre collapses duplicate gene symbols by row-wise mean", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)

  duplicated_df <- data.frame(
    Gene = c("RAC1", "RAC1", "CYFIP1"),
    Sample_01 = c(1, 3, 5),
    Sample_02 = c(2, 4, 6),
    stringsAsFactors = FALSE
  )
  utils::write.csv(duplicated_df, tmp, row.names = FALSE)

  result <- data_pre(tmp)

  expect_equal(nrow(result), 2L)
  expect_equal(unname(result["RAC1", "Sample_01"]), 2)
  expect_equal(unname(result["RAC1", "Sample_02"]), 3)
})

test_that("predict_drs returns the required result columns", {
  path <- example_expression_path()
  input_data <- data_pre(path)
  result <- predict_drs(input_data)

  expect_identical(
    colnames(result),
    c("Sample", "DRS_score", "DRS_group", "n_matched_genes", "n_missing_genes", "missing_genes")
  )
  expect_equal(nrow(result), 6L)
  expect_true(is.numeric(result$DRS_score))
  expect_true(all(result$n_matched_genes == 20L))
  expect_true(all(result$n_missing_genes == 0L))
  expect_true(all(result$missing_genes == ""))
})

test_that("predict_drs errors when model genes are missing", {
  path <- example_expression_path()
  input_data <- data_pre(path)
  missing_input <- input_data[setdiff(rownames(input_data), "RAC1"), , drop = FALSE]

  expect_error(
    predict_drs(missing_input),
    "Missing required model genes: RAC1"
  )
})

test_that("predict_drs errors on single-sample input", {
  path <- example_expression_path()
  input_data <- data_pre(path)[, 1, drop = FALSE]

  expect_error(
    predict_drs(input_data),
    "At least two samples are required"
  )
})

test_that("predict_drs errors on zero-variance model genes", {
  path <- example_expression_path()
  input_data <- data_pre(path)
  input_data["RAC1", ] <- 5

  expect_error(
    predict_drs(input_data),
    "Zero-variance model genes detected: RAC1"
  )
})

test_that("data_pre errors on non-numeric expression values", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)

  non_numeric_df <- data.frame(
    Gene = c("RAC1", "CYFIP1"),
    Sample_01 = c("1.0", "abc"),
    Sample_02 = c("2.0", "3.0"),
    stringsAsFactors = FALSE
  )
  utils::write.csv(non_numeric_df, tmp, row.names = FALSE)

  expect_error(
    data_pre(tmp),
    "Non-numeric expression values detected: abc"
  )
})
