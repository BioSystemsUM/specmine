\name{filter.feature.selection}
\alias{filter.feature.selection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Perform selection by filter
}
\description{
Perform selection by filter using univariate filters, from caret's package.
}
\usage{
filter.feature.selection(datamat, samples.class, functions = rfSBF, method = "cv", repeats = 5)
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
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
caret's sbf result.
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
  ## Example of selection by filter
  rfe.result = filter.feature.selection(metabolomics.dataset$data, "class", functions = rfSBF, method = "repeatedcv", repeats = 10)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~sbf }
\keyword{ ~filters }% __ONLY ONE__ keyword per line
