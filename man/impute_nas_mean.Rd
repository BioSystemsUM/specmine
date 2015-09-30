\name{impute_nas_mean}
\alias{impute_nas_mean}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Impute missing values with mean
}
\description{
  Impute missing values with mean
}
\usage{
impute_nas_mean(dataset)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
      list representing the dataset from a metabolomics experiment.
}
}

\value{
  Returns the dataset with no missing values.
}

\examples{
\dontrun{
  ## Example of NA imputation with mean
  dataset = impute_nas_mean(dataset)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ missing }
\keyword{ values }% __ONLY ONE__ keyword per line
