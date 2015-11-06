\name{metadata_as_variables}
\alias{metadata_as_variables}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Metadata as variables
}
\description{
  Use one or more metadata variables as variables.
}
\usage{
  metadata_as_variables(dataset, metadata.vars, by.index = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{metadata.vars}{
    name or index of the metadata variables that are going to be used as data  variables.
}
  \item{by.index}{
    boolean value indicating if the metadata variables are indexes or names
}
}

\value{
  Returns the dataset with the metadata variables removed from the metadata and added on the data matrix.
}

\examples{
  ## Example of using a metadata variable as data variable
  data(propolis)
  propolis = metadata_as_variables(propolis, "seasons", by.index = FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ metadata }
\keyword{ variable }% __ONLY ONE__ keyword per line
