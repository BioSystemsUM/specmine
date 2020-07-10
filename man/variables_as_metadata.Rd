\name{variables_as_metadata}
\alias{variables_as_metadata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Variables as metadata
}
\description{
  Use one or more data variables as metadata variables.
}
\usage{
  variables_as_metadata(dataset, variables, by.index = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{variables}{
    name or index of the variables that are going to be used as metadata variables.
}
  \item{by.index}{
    boolean value indicating if the variables are indexes or names
}
}

\value{
  Returns the dataset with the variables removed from the data and added on the metadata.
}

\examples{
  ## Example of using a variable as metadata variable
  data(cachexia)
  dataset = variables_as_metadata(cachexia, "Creatine", by.index = FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ metadata }
\keyword{ variable }% __ONLY ONE__ keyword per line
