# DisulfRisk

## Description

`DisulfRisk` is an R package for calculating a 20-gene disulfidptosis-related weighted DRS score from expression matrices.

The package follows a simple user workflow:

1. Read and standardize an expression matrix with `data_pre()`.
2. Load the bundled 20-gene coefficient model with `load_drs_model()`.
3. Calculate sample-level DRS scores with `predict_drs()` or the convenience wrapper `drs_pre()`.
4. Visualize the score distribution with `plot_drs_distribution()`.

`DisulfRisk` is designed for matrices where rows are genes and columns are samples. Duplicate gene symbols are merged by row-wise mean before scoring.

## Model

The DRS score is calculated as:

`DRS = sum(zscore(expression_gene) * coefficient_gene)`

For each model gene, `DisulfRisk` performs a gene-wise z-score across the input sample cohort by default. The bundled model contains the following 20 genes:

`RAC1, CYFIP1, WASF2, BRK1, ABI1, ACTR2, ARPC2, FLNA, FLNB, ACTN4, IQGAP1, MYH10, OSTC, ACSL4, DDO, PDLIM1, DSTN, VCL, CD44, YTHDC1`

## Installation

From the package source directory:

```r
install.packages(".", repos = NULL, type = "source")
```

## Example

```r
library(DisulfRisk)

path <- system.file("extdata", "example_expression.csv", package = "DisulfRisk", mustWork = TRUE)
input_data <- data_pre(path)
result <- drs_pre(input_data)
head(result)
```

## Output

`predict_drs()` returns a data frame with these columns:

- `Sample`
- `DRS_score`
- `DRS_group`
- `n_matched_genes`
- `n_missing_genes`
- `missing_genes`

`DRS_group` is assigned using a default cutoff of `0`: scores greater than or equal to `0` are labeled `High`, and lower scores are labeled `Low`.

## Defensive checks

`DisulfRisk` will stop with an informative error when:

- required model genes are missing
- the expression matrix is not numeric
- only one sample is supplied
- any model gene has zero variance across samples

## Included files

- `inst/extdata/drs_coefficients.csv`: bundled 20-gene coefficient table
- `inst/extdata/example_expression.csv`: example expression matrix for package examples and tests
