\name{x_values_to_indexes}
\alias{x_values_to_indexes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Get x-values indexes
}
\description{
  Returns the indexes corresponding to a vector of x-values (only to numerical values - spectra)
}
\usage{
x_values_to_indexes(dataset, x.values)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
list representing the dataset from a metabolomics experiment.
}
  \item{x.values}{
  vector of x-values.
}
}

\value{
  Returns a numeric vector with the indexes of the x-values.
}

\examples{
  ## Example of getting the indexes of the x-values 
  data(propolis)
  indexes = x_values_to_indexes(propolis, c(1.3, 3.51, 6.03))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ indexes }
\keyword{ xvalues }% __ONLY ONE__ keyword per line
