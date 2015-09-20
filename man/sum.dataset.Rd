\name{sum.dataset}
\alias{sum.dataset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Dataset summary
}
\description{
  Returns a summary with its main features.
}
\usage{
sum.dataset(dataset, stats = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{stats}{
    if TRUE prints some global statistics of the data values.
}
}
\value{
  Returns the summary of the dataset containing:
    \itemize{
      \item Description
      \item Type of data
      \item Number of samples
      \item Number of datapoints
      \item Number of metadata variables if metadata not null
      \item Labels of x axis values and data points if labels not null
    }
  If the parameter 'stats' is TRUE then returns also:
    \itemize{
      \item Number of missing values in data
      \item Mean of data values
      \item Median of data values
      \item Standard deviation
      \item Range of values
      \item Quantiles
    }
}

\examples{
\dontrun{
  ## Example of summarizing a dataset
  sum.dataset(dataset, stats = T)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~dataset }
\keyword{ ~summary }% __ONLY ONE__ keyword per line
