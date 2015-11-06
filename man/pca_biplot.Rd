\name{pca_biplot}
\alias{pca_biplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  PCA biplot
}
\description{
  Shows a PCA biplot.
}
\usage{
pca_biplot(dataset, pca.result, cex = 0.8, legend.cex = 0.8, 
x.colors = 1, inset = c(0, 0), legend.place = "topright", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{pca.result}{
    prcomp object with the PCA results.
}
  \item{cex}{
    cex value.
}
  \item{legend.cex}{
     cex value of the legend.
}
  \item{x.colors}{
	colors of a metadata's variable.
}
  \item{inset}{
    inset parameter of legend function.
}
  \item{legend.place}{
    legend place.
}
  \item{\dots}{
    additional parameters passed to biplot function.
}
}

\examples{
  ## Example of a PCA biplot
  data(cachexia)
  pca.result = pca_analysis_dataset(cachexia)
  pca_biplot(cachexia, pca.result, cex = 0.8)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ pca }
\keyword{ biplot }% __ONLY ONE__ keyword per line
