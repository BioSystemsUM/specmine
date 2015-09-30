\name{pca_scoresplot3D}
\alias{pca_scoresplot3D}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  3D PCA scores plot
}
\description{
    Shows a 3D PCA scores plot of three principal componets.
}
\usage{
pca_scoresplot3D(dataset, pca.result, column.class, 
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
    vector with the principal components to be plotted.
}
}

\examples{
\dontrun{
  ## Example of a 3D PCA scores plot
  pca_scoresplot3D(dataset, pca.result, "metadata.var1", pcas = c(1,2,3))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ scoresplot }
\keyword{ pca }% __ONLY ONE__ keyword per line
