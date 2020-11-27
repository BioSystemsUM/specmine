\name{get_metadata}
\alias{get_metadata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Get metadata
}
\description{
  Get the metadata from the dataset
}
\usage{
get_metadata(dataset)
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
  ## Example of getting the metadata
  library(specmine.datasets)
  data(cachexia)
  cachexia.mt = get_metadata(cachexia)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ metadata }
\keyword{ dataset }% __ONLY ONE__ keyword per line
