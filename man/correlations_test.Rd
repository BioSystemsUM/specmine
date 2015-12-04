\name{correlations_test}
\alias{correlations_test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Correlations test
}
\description{
 Performs correlations test to the whole dataset.
}
\usage{
correlations_test(dataset, method = "pearson", by.var = T, 
alternative = "two.sided")
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
  \item{alternative}{
    alternative argument from cor.test of stats package. Can be "two.sided", "less" or "greater".
}
}

\value{
 Returns a matrix with the correlation values and the p-values
}

\examples{
\dontrun{
  ## Example of correlations test of variables (computationally heavy)
  data(cachexia)
  corr.result = correlations_test(cachexia, 
		method = "pearson", by.var = FALSE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ correlation }
\keyword{ test }% __ONLY ONE__ keyword per line
