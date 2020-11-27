\name{count_missing_values}
\alias{count_missing_values}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Count missing values
}
\description{
  Counts the missing values on the dataset.
}
\usage{
count_missing_values(dataset)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
}

\value{
  Returns the number of missing values on the dataset.
}
\examples{
  ## Example of counting the missing values
  library(specmine.datasets)
  data(cachexia)
  count_missing_values(cachexia)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ missing }
\keyword{ values }% __ONLY ONE__ keyword per line
