\name{subset_metadata}
\alias{subset_metadata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Subset metadata
}
\description{
 Subsets the metadata according to the specified metadata's variables.
}
\usage{
subset_metadata(dataset, variables)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{variables}{
   metadata's variables.
}
}

\value{
 Returns the dataset with the metadata subsetted.
}

\examples{
  ## Example of subsetting samples
  library(specmine.datasets)
  data(propolis)
  subset = subset_metadata(propolis, c("seasons"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ subset }
\keyword{ metadata }% __ONLY ONE__ keyword per line
