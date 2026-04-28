\name{predict_drs}
\alias{predict_drs}
\title{Predict weighted DRS scores for each sample}
\usage{
predict_drs(input_data, model = NULL, group_cutoff = 0)
}
\arguments{
  \item{input_data}{A prepared numeric matrix, a data frame, or a path to an
  expression matrix file.}

  \item{model}{Optional coefficient data frame returned by
  \code{load_drs_model()}. When \code{NULL}, the bundled model is used.}

  \item{group_cutoff}{Numeric cutoff used to define \code{DRS_group}. Scores
  greater than or equal to the cutoff are labeled \code{"High"}; lower scores
  are labeled \code{"Low"}. Defaults to \code{0}.}
}
\value{
A data frame with columns \code{Sample}, \code{DRS_score}, \code{DRS_group},
\code{n_matched_genes}, \code{n_missing_genes}, and \code{missing_genes}.
}
\description{
\code{predict_drs()} applies the bundled 20-gene coefficient model to an
expression matrix. Each model gene is z-scored across the input cohort
before the weighted sum is calculated.
}
\examples{
path <- system.file("extdata", "example_expression.csv", package = "DisulfRisk",
  mustWork = TRUE)
input_data <- data_pre(path)
result <- predict_drs(input_data)
head(result)
}
