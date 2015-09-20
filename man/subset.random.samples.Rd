\name{subset.random.samples}
\alias{subset.random.samples}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Subset random samples
}
\description{
  Gets a subset of random samples from the dataset.
}
\usage{
subset.random.samples(dataset, nsamples)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{nsamples}{
    integer representing the number of samples that we want to get.
}
}

\value{
  Returns the dataset with a number of random samples.
}

\examples{
\dontrun{
  ## Example of subsetting random samples
  subset = subset.random.samples(dataset, 10)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~random }
\keyword{ ~subset }% __ONLY ONE__ keyword per line
