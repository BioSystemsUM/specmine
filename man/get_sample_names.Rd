\name{get_sample_names}
\alias{get_sample_names}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Get sample names
}
\description{
  Get the sample names from the dataset.
}
\usage{
get_sample_names(dataset)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
}

\value{
  Returns a vector with the sample names.
}

\examples{
  ## Example of getting the sample names
  library(specmine.datasets)
  data(cachexia)
  sample.names = get_sample_names(cachexia)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ sample }
\keyword{ dataset }% __ONLY ONE__ keyword per line
