\name{feature_selection}
\alias{feature_selection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Perform feature selection
}
\description{
  Perform feature selection on the dataset.
}
\usage{
feature_selection(dataset, column.class, method = "rfe", 
functions, validation = "cv", repeats = 5, number = 10, 
subsets = 2^(2:4))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{column.class}{
    string or index indicating what metadata to use.
}
  \item{method}{
    method used for feature selection. Possible values are "rfe" (recursive feature elimination) and "filter" (Selection by filter - sbf) from caret's package.
}
  \item{functions}{
    a list of functions for model fitting, prediction and variable importance/filtering.
}
  \item{validation}{
    the external resampling method: boot, cv, LOOCV or LGOCV (for repeated training/test splits.
}
  \item{repeats}{
    for repeated k-fold cross-validation only: the number of complete sets of folds to compute.
}
  \item{number}{
    either the number of folds or number of resampling iterations.
}
  \item{subsets}{
    a numeric vector of integers corresponding to the number of features that should be retained (rfe only).
}
}

\value{
  caret's result from rfe or sbf.
}

\examples{
\donttest{
  ## Example of feature selection using rfe and sbf
  library(caret)
  rfe.result = feature_selection(cachexia, "Muscle.loss", 
               method="rfe", functions = caret::rfFuncs, 
               validation = "cv", number = 3, 
               subsets = 2^(1:6))
  sbf.result = feature_selection(cachexia, "Muscle.loss", 
	       method="filter", functions = caret::rfSBF, 
	       validation = "cv")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ featureselection }
\keyword{ filters }% __ONLY ONE__ keyword per line
\keyword{ wrappers }
