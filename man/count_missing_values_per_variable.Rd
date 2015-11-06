\name{count_missing_values_per_variable}
\alias{count_missing_values_per_variable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Count missing values per variable
}
\description{
  Counts the missing values on each variable of the dataset.
}
\usage{
count_missing_values_per_variable(dataset, remove.zero = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{remove.zero}{
    boolean value indicating if the results of variables with no missing value are removed.
}
}

\value{
  Returns a vector with the number of missing values on each sample.
}
\examples{
  ## Example of counting the missing values on each variable
  data(cachexia)
  cachexia$data[10,10] = NA
  count_missing_values_per_variable(cachexia)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ missing }
\keyword{ values }% __ONLY ONE__ keyword per line
