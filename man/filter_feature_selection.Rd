\name{filter_feature_selection}
\alias{filter_feature_selection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Perform selection by filter
}
\description{
Perform selection by filter using univariate filters, from caret's package.
}
\usage{
filter_feature_selection(datamat, samples.class, 
functions = caret::rfSBF, method = "cv", repeats = 5)
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
    a list of functions for model fitting, prediction and variable filtering.
}
  \item{method}{
    the external resampling method: boot, cv, LOOCV or LGOCV (for repeated training/test splits.
}
  \item{repeats}{
    for repeated k-fold cross-validation only: the number of complete sets of folds to compute.
}
}

\value{
A caret's sbf object with the result of selection by filter.
}

\examples{
\donttest{
  ## Example of selection by filter
  library(caret)
  library(specmine.datasets)
  data(cachexia)
  rfe.result = filter_feature_selection(cachexia$data, 
	       cachexia$metadata$Muscle.loss, functions = caret::rfSBF, 
	       method = "cv")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ sbf }
\keyword{ filters }% __ONLY ONE__ keyword per line
