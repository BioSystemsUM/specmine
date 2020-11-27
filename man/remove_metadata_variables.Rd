\name{remove_metadata_variables}
\alias{remove_metadata_variables}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Remove metadata's variables
}
\description{
 Remove metadata's variables from the dataset
}
\usage{
remove_metadata_variables(dataset, variables.to.remove)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{variables.to.remove}{
    vector with the metadata's variables to remove.
}
}

\value{
  Returns the dataset with the selected metadata's variables removed.
}


\examples{
  ## Example of removing metadata's variables
  library(specmine.datasets)
  data(propolis)
  dataset = remove_metadata_variables(propolis, c("seasons"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ remove }
\keyword{ metadata }% __ONLY ONE__ keyword per line
