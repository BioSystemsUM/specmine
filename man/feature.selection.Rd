\name{feature.selection}
\alias{feature.selection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Perform feature Selection
}
\description{
  Perform feature selection on the dataset.
}
\usage{
feature.selection(dataset, column.class, method = "rfe", functions, validation = "cv", repeats = 5, number = 10, subsets = 2^(2:4))
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
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  caret's result from rfe or sbf.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
\dontrun{
  ## Example of feature selection using rfe and sbf
  rfe.result = feature.selection(metabolomics.dataset, "class", method="rfe",functions = rfFuncs, validation = "repeatedcv", repeats = 10, number = 5, subsets = 2^(1:6))
  sbf.result = feature.selection(metabolomics.dataset, "class", method="filter", functions = rfSBF, validation = "repeatedcv", repeats = 5)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~featureselection }
\keyword{ ~filters }% __ONLY ONE__ keyword per line
\keyword{ ~wrappers }