\name{recursive.feature.elimination}
\alias{recursive.feature.elimination}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Perform recursive feature elimination
}
\description{
Perform recursive feature elimination on the dataset using caret's package.
}
\usage{
recursive.feature.elimination(datamat, samples.class, functions = rfFuncs, method = "cv", repeats = 5, number = 10, subsets = 2^(2:4))
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
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A caret's rfe object with the result of recursive feature selection.
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
  ## Example of recursive feature elimination
  rfe.result = recursive.feature.elimination(metabolomics.dataset$data, "class", functions = rfFuncs, method = "repeatedcv", repeats = 10, number = 5, subsets = 2^(1:6))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~rfe }
\keyword{ ~wrappers }% __ONLY ONE__ keyword per line
