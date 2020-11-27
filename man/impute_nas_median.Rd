\name{impute_nas_median}
\alias{impute_nas_median}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Impute missing values with median
}
\description{
  Impute missing values with median
}
\usage{
impute_nas_median(dataset)
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
  ## Example of NA imputation with median
  library(specmine.datasets)
  data(propolis)
  propolis = impute_nas_median(propolis)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ missing }
\keyword{ values }% __ONLY ONE__ keyword per line
