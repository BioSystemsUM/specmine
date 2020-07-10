\name{merge_data_metadata}
\alias{merge_data_metadata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Merge data and metadata
}
\description{
  Merges the data and metadata from the dataset into a single data.frame.
}
\usage{
merge_data_metadata(dataset, samples = NULL, 
metadata.vars = NULL, x.values = NULL, by.index = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{samples}{
    vector with indexes or names of the samples to select
}
  \item{metadata.vars}{
    metadata's variables.
}
  \item{x.values}{
    vector with the desired x-values to get from the dataset.
}
  \item{by.index}{
    if TRUE, the values of the x.values argument are indexes.  
}
}

\value{
 Returns a data.frame with the data and metadata from the dataset merged.
}

\examples{
  ## Example of merging data and metadata
  data(cachexia)
  dt.merged = merge_data_metadata(cachexia)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ data }
\keyword{ metadata }% __ONLY ONE__ keyword per line
