\name{get.sample.names}
\alias{get.sample.names}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Get sample names
}
\description{
  Get the sample names from the dataset.
}
\usage{
get.sample.names(dataset)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
}

\value{
  Returns a vector with the sample names.
}

\examples{
\dontrun{
  ## Example of getting the sample names
  sample.names = get.sample.names(dataset)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~sample }
\keyword{ ~dataset }% __ONLY ONE__ keyword per line
