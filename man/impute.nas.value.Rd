\name{impute.nas.value}
\alias{impute.nas.value}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Impute missing values with value replacement
}
\description{
  Impute missing values with value replacement.
}
\usage{
impute.nas.value(dataset, value)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
      list representing the dataset from a metabolomics experiment.
}
  \item{value}{
    value to replace the missing values.
}
}

\value{
  Returns the dataset with no missing values.
}

\examples{
\dontrun{
  ## Example of NA imputation with value replacing
  dataset = impute.nas.value(dataset, 0.0005)
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~missing }
\keyword{ ~values }% __ONLY ONE__ keyword per line
