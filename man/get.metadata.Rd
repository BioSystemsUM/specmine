\name{get.metadata}
\alias{get.metadata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Get metadata
}
\description{
  Get the metadata from the dataset
}
\usage{
get.metadata(dataset)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
}
\value{
  returns a data frame with the metadata.
}

\examples{
\dontrun{
  ## Example of getting the metadata
  metadata = get.metadata(dataset)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~metadata }
\keyword{ ~dataset }% __ONLY ONE__ keyword per line
