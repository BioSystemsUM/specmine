\name{subset_x_values}
\alias{subset_x_values}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Subset x-values
}
\description{
 Gets a subset of specific x-values from the dataset.
}
\usage{
subset_x_values(dataset, variables, by.index = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{variables}{
    vector with the desired variables to get from the dataset.
}
  \item{by.index}{
    if TRUE, the values of the variables argument are indexes.  
}
}

\value{
  Returns the dataset with the selected set of x-values.
}

\examples{
  ## Example of subsetting x-values
  library(specmine.datasets)
  data(cachexia)
  subset = subset_x_values(cachexia, c(1,2,10,20), by.index = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ subset }
\keyword{ xvalues }% __ONLY ONE__ keyword per line
