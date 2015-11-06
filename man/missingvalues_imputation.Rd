\name{missingvalues_imputation}
\alias{missingvalues_imputation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Missing values imputation
}
\description{
  Treats the missing values of a dataset according to a specific method.
}
\usage{
missingvalues_imputation(dataset, method = "value", 
value = 5e-04, k = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{method}{
    imputation method. It can be:
      \itemize{
        \item \bold{"value"} - replaces the missing values with a specific value
        \item \bold{"mean"} - replaces the missing values with the mean of the variables' values
        \item \bold{"median"} - replaces the missing values with the median of the variables' values
        \item \bold{"knn"} - replaces the missing values with k nearest neighbor averaging
        \item \bold{"linapprox"} - replaces the missing values with linear approximation
      }
}
  \item{value}{
    the value to replace the missing values if the method is "value".
}
  \item{k}{
    the number of neighbors if the method is "knn".
}
}

\value{
  Returns the dataset with no missing values.
}
\examples{
  ## Example of impute missing values
  data(propolis)
  dataset = missingvalues_imputation(propolis, method = "value", 
	    value = 0.0005)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ missing }
\keyword{ values }% __ONLY ONE__ keyword per line
