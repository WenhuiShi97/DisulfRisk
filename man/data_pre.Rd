\name{data_pre}
\alias{data_pre}
\title{Prepare an expression matrix for DRS scoring}
\usage{
data_pre(x, gene_col = 1, sep = NULL)
}
\arguments{
  \item{x}{A file path to a csv/tsv/txt matrix, or an in-memory matrix/data
  frame.}

  \item{gene_col}{For file or data frame inputs where gene symbols are stored
  in a column, the 1-based index of that column. Defaults to \code{1}.}

  \item{sep}{Optional delimiter override. By default the delimiter is inferred
  from the file extension or file contents.}
}
\value{
A numeric matrix with gene symbols as row names and samples as columns.
}
\description{
\code{data_pre()} reads a gene-by-sample expression matrix from a file or
standardizes an in-memory matrix/data frame. By default, rows are genes
and columns are samples. Duplicate gene symbols are merged by row-wise mean.
}
\examples{
path <- system.file("extdata", "example_expression.csv", package = "DisulfRisk",
  mustWork = TRUE)
input_data <- data_pre(path)
dim(input_data)
}
