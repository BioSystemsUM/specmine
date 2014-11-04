\name{correlations.dataset}
\alias{correlations.dataset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Dataset correlations
}
\description{
  Calculate the correlations of all variables or samples in the dataset.
}
\usage{
correlations.dataset(dataset, method = "pearson", by.var = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{method}{
    correlation method, it can be "pearson", "kendall" or "spearman".
}
  \item{by.var}{
    if TRUE then the correlations of the variables will be calculated, if not then the correlations of the samples will be calculated.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  correlation matrix.
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
  ## Example of correlations of variables
  corr.result = correlations.dataset(metabolomics.dataset, method = "pearson", by.var = T)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~correlation }
\keyword{ ~spearman }% __ONLY ONE__ keyword per line
\keyword{ ~kendall }% __ONLY ONE__ keyword per line
\keyword{ ~pearson }% __ONLY ONE__ keyword per line
