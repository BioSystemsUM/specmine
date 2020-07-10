\name{correlation_test}
\alias{correlation_test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Correlation test of two variables or samples
}
\description{
 Performs correlations test of two variables or samples from the dataset.
}
\usage{
correlation_test(dataset, x, y, method = "pearson", 
alternative = "two.sided", by.var = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{x}{
    first variable or sample.
}
  \item{y}{
   second variable or sample.
}
  \item{method}{
    correlation method, it can be "pearson", "kendall" or "spearman".
}
  \item{alternative}{
    alternative argument from cor.test of stats package. Can be "two.sided", "less" or "greater".
}
  \item{by.var}{
    if TRUE then the correlations of the variables will be calculated, if not then the correlations of the samples will be calculated.
}
}

\value{
 Returns the correlation result from cor.test function of stats package.
}

\examples{
  ## Example of correlations test of variables
  data(cachexia)
  corr.result = correlation_test(cachexia, "Serine", "Creatine", method = "pearson", 
		by.var = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ correlation }
\keyword{ test }% __ONLY ONE__ keyword per line
