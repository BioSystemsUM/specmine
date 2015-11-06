\name{pca_importance}
\alias{pca_importance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  PCA importance
}
\description{
  Gets the importance from the PCs.
}
\usage{
pca_importance(pca.res, pcs = 1:length(pca.res$sdev), sd = T, 
prop = T, cumul = T, min.cum = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pca.res}{
    prcomp object with the PCA results.
}
  \item{pcs}{
    vector with the PCs to get.
}
  \item{sd}{
    boolean value indicating if standard deviation will be returned or not.
}
  \item{prop}{
    boolean value that indicates if the proportion of variance is returned or not.
}
  \item{cumul}{
    boolean value that indicates if the cumulative variance is returned or not.
}
  \item{min.cum}{
    allows to define minimum cumulative \% of variance
}
}

\value{
  Returns the information about the importance of the PCs.
}

\examples{
  ## Example of performing a classical PCA analysis
  data(cachexia)
  pca.result = pca_analysis_dataset(cachexia)
  pca_importance(pca.result, pcs = 1:5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ pca }
\keyword{ importance }% __ONLY ONE__ keyword per line
