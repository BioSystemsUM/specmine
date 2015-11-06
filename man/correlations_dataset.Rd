\name{correlations_dataset}
\alias{correlations_dataset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Dataset correlations
}
\description{
  Calculate the correlations of all variables or samples in the dataset.
}
\usage{
correlations_dataset(dataset, method = "pearson", by.var = T)
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

\value{
  Returns the correlation matrix
}


\examples{
  ## Example of correlations of variables
  data(cachexia)
  corr.result = correlations_dataset(cachexia, 
		method = "pearson", by.var = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ correlation }
\keyword{ spearman }% __ONLY ONE__ keyword per line
\keyword{ kendall }% __ONLY ONE__ keyword per line
\keyword{ pearson }% __ONLY ONE__ keyword per line
