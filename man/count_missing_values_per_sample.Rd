\name{count_missing_values_per_sample}
\alias{count_missing_values_per_sample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Count missing values per sample
}
\description{
  Counts the missing values on each sample of the dataset.
}
\usage{
count_missing_values_per_sample(dataset, remove.zero = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{remove.zero}{
    boolean value indicating if the results of samples with no missing value are removed.
}
}

\value{
  Returns a vector with the number of missing values on each sample.
}
\examples{
  ## Example of counting the missing values on each sample
  library(specmine.datasets)
  data(cachexia)
  cachexia$data[10,10] = NA
  count_missing_values_per_sample(cachexia)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ missing }
\keyword{ values }% __ONLY ONE__ keyword per line
