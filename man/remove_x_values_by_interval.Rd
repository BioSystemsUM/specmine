\name{remove_x_values_by_interval}
\alias{remove_x_values_by_interval}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Remove x-values by interval
}
\description{
  Remove an interval of x-values from the dataset.
}
\usage{
remove_x_values_by_interval(dataset, min.value, max.value)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{min.value}{
    minimum value of the interval.
}
  \item{max.value}{
    maximum value of the interval.
}
}

\value{
 Returns the dataset with the specified interval of x-values removed.
}

\examples{
  ## Example of removing x-values by interval
  data(cassavaPPD)
  cassavaPPD = remove_x_values_by_interval(cassavaPPD, 200, 300)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ remove }
\keyword{ xvalue }% __ONLY ONE__ keyword per line
