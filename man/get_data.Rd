\name{get_data}
\alias{get_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Get data
}
\description{
  Get the data matrix from dataset
}
\usage{
get_data(dataset)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
}

\value{
  Returns the data matrix
}

\examples{
  ## Example of getting the data matrix
  data(cachexia)
  cachexia.dm = get_data(cachexia)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ matrix }
\keyword{ dataset }% __ONLY ONE__ keyword per line
