\name{xvalue_interval_to_indexes}
\alias{xvalue_interval_to_indexes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Get indexes of an interval of x-values
}
\description{
Returns indexes corresponding to an interval of x-values (only to numerical values - spectra)
}
\usage{
xvalue_interval_to_indexes(dataset, min.value, max.value)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
  list representing the dataset from a metabolomics experiment.
}
  \item{min.value}{
mininum x-value of the interval.
}
  \item{max.value}{
maximum x-value of the interval.
}
}

\value{
Returns a numeric vector with the indexes of the x-values interval
}

\examples{
  ## Example of getting the indexes of an interval of x-values
  library(specmine.datasets)
  data(propolis)
  indexes.interval = xvalue_interval_to_indexes(propolis, 2.0, 5.5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ indexes }
\keyword{ interval }% __ONLY ONE__ keyword per line
