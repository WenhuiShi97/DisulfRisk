\name{drs_pre}
\alias{drs_pre}
\title{Convenience wrapper for DRS prediction}
\usage{
drs_pre(input_data, model = NULL, group_cutoff = 0)
}
\arguments{
  \item{input_data}{A prepared numeric matrix, a data frame, or a path to an
  expression matrix file.}

  \item{model}{Optional coefficient data frame returned by
  \code{load_drs_model()}. When \code{NULL}, the bundled model is used.}

  \item{group_cutoff}{Numeric cutoff used to define \code{DRS_group}. Defaults
  to \code{0}.}
}
\value{
A data frame with DRS prediction results.
}
\description{
\code{drs_pre()} is a lightweight wrapper around \code{\link{predict_drs}}
that keeps the user workflow concise in examples and scripts.
}
\examples{
path <- system.file("extdata", "example_expression.csv", package = "DisulfRisk",
  mustWork = TRUE)
input_data <- data_pre(path)
result <- drs_pre(input_data)
head(result)
}
