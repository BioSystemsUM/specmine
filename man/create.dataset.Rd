\name{create.dataset}
\alias{create.dataset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Create dataset
}
\description{
  Create a dataset from existing objects
}
\usage{
create.dataset(datamatrix, type = "undefined", metadata = NULL, description = "", sample.names = NULL, x.axis.values = NULL, label.x = NULL, label.values = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datamatrix}{
  matrix with numerical data: rows are assumed to be variables and columns assumed to be samples.
}
  \item{type}{
  type of data: string that can be one of the following:
  \itemize{
    \item nmr-spectra
    \item nmr-peaks
    \item ir-spectra
    \item uvv-spectra
    \item concentrations
    \item undefined
  }
  \item{metadata}{
  data frame with the dataset's metadata: columns represent each metadata variable and rows represent the value of the metadata for the sample.
}
  \item{description}{
  string with a short description of the dataset.
}
  \item{sample.names}{
  vector with sample names, if NULL then the column names of datamatrix or sequential numbers will be used.
}
  \item{x.axis.values}{
  vector with the x axis values, if NULL then the row names of datamatrix or sequential numbers will be used.
}
  \item{label.x}{
  x axis label.
}
  \item{label.values}{
  values label.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  list representing the dataset:
    \item{data }{matrix with the data}
    \item{type }{type of the data}
    \item{description }{short description of the dataset}
    \item{metadata }{data frame with the metadata variables}
    \item{labels }{list with labels of x axis and values}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
  \dontrun{
  ## Example of creating a dataset
  dataset = create.dataset(datamatrix, type = "nmr-spectra", metadata = dataset.metadata, description = "short description", sample.names = NULL, x.axis.values = NULL, label.x = "ppm", label.values = "intensity")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~dataset }
