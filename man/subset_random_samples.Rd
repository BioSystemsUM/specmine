\name{subset_random_samples}
\alias{subset_random_samples}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Subset random samples
}
\description{
  Gets a subset of random samples from the dataset.
}
\usage{
subset_random_samples(dataset, nsamples)
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
  ## Example of subsetting random samples
  data(cachexia)
  subset = subset_random_samples(cachexia, 15)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ random }
\keyword{ subset }% __ONLY ONE__ keyword per line
