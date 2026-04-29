\name{plot_drs_distribution}
\alias{plot_drs_distribution}
\title{Plot the DRS score distribution}
\usage{
plot_drs_distribution(
  x,
  model = NULL,
  group_cutoff = 0,
  main = "DRS score distribution",
  ylab = "DRS score",
  ylim = NULL,
  las = 2,
  ...
)
}
\arguments{
  \item{x}{A prediction data frame returned by \code{predict_drs()}, or any
  input accepted by \code{drs_pre()}.}

  \item{model}{Optional coefficient data frame returned by
  \code{load_drs_model()} when \code{x} is not already a prediction result.}

  \item{group_cutoff}{Numeric cutoff used to define \code{DRS_group} when
  \code{x} is not already a prediction result. Defaults to \code{0}.}

  \item{main}{Plot title.}

  \item{ylab}{Y-axis label.}

  \item{ylim}{Optional numeric vector of length 2 giving the y-axis range.}

  \item{las}{Sample label orientation. The default value \code{2} draws labels
  vertically.}

  \item{...}{Additional arguments passed to \code{\link[graphics]{barplot}}.}
}
\value{
Invisibly returns the plotted prediction data frame.
}
\description{
\code{plot_drs_distribution()} draws a simple base-R bar plot of sample-level
DRS scores ordered from low to high. The function can accept a prediction
result from \code{\link{predict_drs}} or any input supported by
\code{\link{drs_pre}}.
}
\examples{
path <- system.file("extdata", "example_expression.csv", package = "DisulfRisk",
  mustWork = TRUE)
input_data <- data_pre(path)
result <- drs_pre(input_data)
\donttest{
  plot_drs_distribution(result)
}
}
