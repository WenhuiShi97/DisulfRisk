\name{load_drs_model}
\alias{load_drs_model}
\title{Load and validate the bundled DRS coefficient model}
\usage{
load_drs_model(file = NULL)
}
\arguments{
  \item{file}{Optional path to a coefficient csv file. When \code{NULL}, the
  bundled \code{drs_coefficients.csv} file is used.}
}
\value{
A data frame with columns \code{Gene} and \code{Coefficient}.
}
\description{
\code{load_drs_model()} reads \code{drs_coefficients.csv} from
\code{inst/extdata/} and validates that the coefficient table contains the
required \code{Gene} and \code{Coefficient} columns for the complete 20-gene
DRS model.
}
\examples{
model <- load_drs_model()
head(model)
}
