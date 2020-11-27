\name{get_data_as_df}
\alias{get_data_as_df}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Get data as data frame
}
\description{
  Get the data matrix from the dataset as a data frame.
}
\usage{
get_data_as_df(dataset)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
}

\value{
  Returns the data matrix from the dataset as a data.frame object.
}

\examples{
  ## Example of getting the data matrix as data frame
  library(specmine.datasets)
  data(cachexia)
  cachexia.dt = get_data_as_df(cachexia)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ matrix }
\keyword{ dataframe }% __ONLY ONE__ keyword per line
