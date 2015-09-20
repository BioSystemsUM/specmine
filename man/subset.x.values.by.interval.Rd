\name{subset.x.values.by.interval}
\alias{subset.x.values.by.interval}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Subset x-values by interval
}
\description{
 Gets a subset of a specific interval of x-values.
}
\usage{
subset.x.values.by.interval(dataset, min.value, max.value)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{min.value}{
    the minimum value of the interval.
}
  \item{max.value}{
    the maximum value of the interval.
}
}

\value{
  Returns the dataset with the selected interval of x-values.
}

\examples{
\dontrun{
  ## Example of subsetting x-values by an interval
  subset = subset.x.values.by.interval(dataset, 200, 800)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~subset }
\keyword{ ~interval }% __ONLY ONE__ keyword per line
