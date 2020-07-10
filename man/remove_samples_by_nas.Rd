\name{remove_samples_by_nas}
\alias{remove_samples_by_nas}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Remove samples by NAs
}
\description{
  Remove samples from the dataset by the number of NAs
}
\usage{
remove_samples_by_nas(dataset, max.nas = 0, by.percent = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{max.nas}{
    number of NAs (or percentage) to which samples with more NAs will be removed.
}
  \item{by.percent}{
    if TRUE the max.nas argument is a percentage.
}
}

\value{
  Returns the dataset with the samples with more NAs than the limit removed.
}

\examples{
  ## Example of removing samples by NAs
  data(propolis)
  propolis = remove_samples_by_nas(propolis, 40, by.percent = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ remove }
\keyword{ missingvalue }% __ONLY ONE__ keyword per line
