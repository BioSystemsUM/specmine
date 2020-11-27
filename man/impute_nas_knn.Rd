\name{impute_nas_knn}
\alias{impute_nas_knn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Impute missing values with KNN
}
\description{
  Impute missing values with KNN
}
\usage{
impute_nas_knn(dataset, k = 10, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
      list representing the dataset from a metabolomics experiment.
}
  \item{k}{
    the number of nearest neighbors.
}
  \item{\dots}{
    additional values to impute.knn function.
}
}

\value{
   Returns the dataset with no missing values.
}

\examples{
  ## Example of NA imputation with knn
  library(specmine.datasets)
  data(propolis)
  dataset = impute_nas_knn(propolis, k=10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ missing }
\keyword{ values }% __ONLY ONE__ keyword per line
