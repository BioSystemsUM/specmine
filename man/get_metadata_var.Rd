\name{get_metadata_var}
\alias{get_metadata_var}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Get metadata variable
}
\description{
  Get the values of a metadata variable from the dataset.
}
\usage{
get_metadata_var(dataset, var)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{var}{
    index or name of the metadata variable.
}
}

\value{
  Returns a vector with the values of the metadata variable. 
}

\examples{
  ## Example of getting a metadata variable
  library(specmine.datasets)
  data(cachexia)
  metadata.variable = get_metadata_var(cachexia, "Muscle.loss")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ metadata }
\keyword{ variable }% __ONLY ONE__ keyword per line
