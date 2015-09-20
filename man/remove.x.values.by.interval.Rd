\name{remove.x.values.by.interval}
\alias{remove.x.values.by.interval}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Remove x-values by interval
}
\description{
  Remove an interval of x-values from the dataset.
}
\usage{
remove.x.values.by.interval(dataset, min.value, max.value)
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
\dontrun{
  ## Example of removing x-values by interval
  dataset = remove.x.values.by.interval(dataset, 200, 300)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~remove }
\keyword{ ~xvalue }% __ONLY ONE__ keyword per line
