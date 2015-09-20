\name{remove.samples}
\alias{remove.samples}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Remove samples
}
\description{
 Remove samples from the dataset.
}
\usage{
remove.samples(dataset, samples.to.remove, rebuild.factors = T)
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
\dontrun{
  ## Example of removing samples
  dataset = remove.samples(dataset, c("sample1","sample2"))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~remove }
\keyword{ ~sample }% __ONLY ONE__ keyword per line
