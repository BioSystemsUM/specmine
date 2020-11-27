\name{get_data_values}
\alias{get_data_values}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get data values
}
\description{
Gets the values of all samples in the dataset given a set of x axis names or indexes.
}
\usage{
get_data_values(dataset, x.axis.val, by.index = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
list representing the dataset from a metabolomics experiment.
}
  \item{x.axis.val}{
vector with the values of the x axis (could be names or indexes).
}
  \item{by.index}{
boolean value indicating if the x.axis.val is a vector of indexes or not.
}
}

\value{
Returns a matrix with the values of all samples in the specified x axis.
}

\examples{
  ## Example of getting a metadata value
  library(specmine.datasets)
  data(cachexia)
  data.values = get_data_values(cachexia, c("Creatine","Serine","Lactate"), 
		by.index = FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ data }
\keyword{ values }% __ONLY ONE__ keyword per line
