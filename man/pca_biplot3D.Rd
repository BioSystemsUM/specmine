\name{pca_biplot3D}
\alias{pca_biplot3D}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  3D PCA biplot (interactive)
}
\description{
  Shows a interactive 3D PCA biplot.
}
\usage{
pca_biplot3D(dataset, pca.result, column.class = NULL,
pcas = c(1, 2, 3))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
        list representing the dataset from a metabolomics experiment.
}
  \item{pca.result}{
    prcomp object with the PCA results.
}
  \item{column.class}{
    metadata's variable.
}
  \item{pcas}{
    the three principal components.
}
}

\examples{
\donttest{
  ### Example of a 3D PCA biplot
  library(specmine.datasets)
  data(cachexia)
  pca.result = pca_analysis_dataset(cachexia)
  pca_biplot3D(cachexia, pca.result, "Muscle.loss", pcas = c(1,2,3))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ pca }
\keyword{ biplot }% __ONLY ONE__ keyword per line
