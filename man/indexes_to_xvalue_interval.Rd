\name{indexes_to_xvalue_interval}
\alias{indexes_to_xvalue_interval}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Get the x-values of a vector of indexes
}
\description{
  Returns x-values corresponding to a vector of indexes (only to numerical values - spectra)
}
\usage{
indexes_to_xvalue_interval(dataset, indexes)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
  list representing the dataset from a metabolomics experiment.
}
  \item{indexes}{
  numeric vector containing the indexes.
}
}

\value{
Returns a numeric vector with the interval of x-values from the indexes vector
}

\examples{
  ## Example of getting the interval of x-values from indexes
  data(cassavaPPD)
  xvalue.interval = indexes_to_xvalue_interval(cassavaPPD, c(10,50))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ xvalues }
\keyword{ indexes }% __ONLY ONE__ keyword per line
