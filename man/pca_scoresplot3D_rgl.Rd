\name{pca_scoresplot3D_rgl}
\alias{pca_scoresplot3D_rgl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  3D PCA scores plot (interactive)
}
\description{
  Shows a interactive 3D PCA scores plot of three principal components.
}
\usage{
pca_scoresplot3D_rgl(dataset, pca.result, column.class, 
pcas = c(1, 2, 3), size = 1, labels = FALSE)
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
  \item{size}{
    parameter of plot3d from rgl package.
}
  \item{labels}{
    boolean value indicating if the samples' labels will be shown.
}
}


\examples{
  ## Example of a 3D PCA scores plot
  data(cachexia)
  pca.result = pca_analysis_dataset(cachexia)
  pca_scoresplot3D_rgl(cachexia, pca.result, "Muscle.loss", 
  		       pcas = c(1,2,3), labels = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ scoresplot }
\keyword{ pca }% __ONLY ONE__ keyword per line
