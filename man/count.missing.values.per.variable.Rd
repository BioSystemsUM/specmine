\name{count.missing.values.per.variable}
\alias{count.missing.values.per.variable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Count missing values per variable
}
\description{
  Counts the missing values on each variable of the dataset.
}
\usage{
count.missing.values.per.variable(dataset, remove.zero = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{remove.zero}{
    boolean value indicating if the results of variables with no missing value are removed.
}
}

\value{
  Returns a vector with the number of missing values on each sample.
}
\examples{
\dontrun{
  ## Example of counting the missing values on each variable
  count.missing.values.per.variable(dataset)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~missing }
\keyword{ ~values }% __ONLY ONE__ keyword per line
