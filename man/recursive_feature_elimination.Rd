\name{recursive_feature_elimination}
\alias{recursive_feature_elimination}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Perform recursive feature elimination
}
\description{
Perform recursive feature elimination on the dataset using caret's package.
}
\usage{
recursive_feature_elimination(datamat, samples.class, 
functions = caret::rfFuncs, method = "cv", repeats = 5, 
number = 10, subsets = 2^(2:4))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datamat}{
    data matrix from dataset.
}
  \item{samples.class}{
    string or index indicating what metadata to use.
}
  \item{functions}{
    a list of functions for model fitting, prediction and variable importance.
}
  \item{method}{
    the external resampling method: boot, cv, LOOCV or LGOCV (for repeated training/test splits.
}
  \item{repeats}{
    for repeated k-fold cross-validation only: the number of complete sets of folds to compute.
}
  \item{number}{
    either the number of folds or number of resampling iterations.
}
  \item{subsets}{
    a numeric vector of integers corresponding to the number of features that should be retained.
}
}

\value{
A caret's rfe object with the result of recursive feature selection.
}

\examples{
\donttest{
  ## Example of recursive feature elimination
  library(specmine.datasets)
  data(cachexia)
  library(caret)
  rfe.result = recursive_feature_elimination(cachexia$data, 
	       cachexia$metadata$Muscle.loss, functions = caret::rfFuncs, 
	       method = "cv", number = 3, subsets = 2^(1:6))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ rfe }
\keyword{ wrappers }% __ONLY ONE__ keyword per line
