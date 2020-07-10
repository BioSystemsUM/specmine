\name{subset_samples}
\alias{subset_samples}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Subset samples
}
\description{
  Gets a subset of specific samples from the dataset.
}
\usage{
subset_samples(dataset, samples, rebuild.factors = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{samples}{
    vector with indexes or names of the samples to select
}
  \item{rebuild.factors}{
    if TRUE the metadata factors are rebuilded.
}
}

\value{
  Returns the dataset with the selected set of samples.
}

\examples{
  ## Example of subsetting samples
  data(cachexia)
  subset = subset_samples(cachexia, c("PIF_178","PIF_132"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ subset }
\keyword{ sample }% __ONLY ONE__ keyword per line
