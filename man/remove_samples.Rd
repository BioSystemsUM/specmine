\name{remove_samples}
\alias{remove_samples}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Remove samples
}
\description{
 Remove samples from the dataset.
}
\usage{
remove_samples(dataset, samples.to.remove, rebuild.factors = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{samples.to.remove}{
    vector with the sample's names to remove.
}
  \item{rebuild.factors}{
    if TRUE, rebuilds the factors from metadata.
}
}

\value{
 Returns the dataset with the specified samples removed. 
}

\examples{
  ## Example of removing samples
  data(cachexia)
  cachexia = remove_samples(cachexia, c("PIF_178","PIF_090"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ remove }
\keyword{ sample }% __ONLY ONE__ keyword per line
