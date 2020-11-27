\name{get_type}
\alias{get_type}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Get type of data
}
\description{
  Get the type of the data from the dataset
}
\usage{
get_type(dataset)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
}

\value{
  Returns a string with the type of the data.
}

\examples{
  ## Example of getting the type of the data
  library(specmine.datasets)
  data(cachexia)
  type = get_type(cachexia)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ type }
\keyword{ dataset }% __ONLY ONE__ keyword per line
