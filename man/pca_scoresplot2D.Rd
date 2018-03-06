\name{pca_scoresplot2D}
\alias{pca_scoresplot2D}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  2D PCA scores plot
}
\description{
  Shows a 2D PCA scores plot of two principal componets.
}
\usage{
pca_scoresplot2D(dataset, pca.result, column.class, 
pcas = c(1, 2), labels = FALSE, ellipses = FALSE, bw=FALSE,
pallette = 2, leg.pos = "right", xlim = NULL, ylim = NULL)
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
    vector of two elements with the PCs that will be plotted.
}
  \item{labels}{
    boolean value indicating if the sample's labels will be displayed.
}
  \item{ellipses}{
    boolean value that indicates if an ellipse will be drawn on each group of the metadata's variable. Ellipses will not be drawn if bw=TRUE.
}
  \item{bw}{
    if TRUE, it will be displayed a black and white plot. It defaults to FALSE.
}
  \item{pallette}{
    parameter of scale_colour_brewer from ggplot2.
}
  \item{leg.pos}{
    position of the legend.
}
  \item{xlim}{
    vector with two numeric values indicating the limits of the x axis.
}
  \item{ylim}{
    vector with two numeric values indicating the limits of the y axis.
}
}

\examples{
  ## Example of a 2D PCA scores plot
  data(cachexia)
  pca.result = pca_analysis_dataset(cachexia)
  pca_scoresplot2D(cachexia, pca.result, "Muscle.loss", pcas = c(1,2), 
    		   ellipses = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ scoresplot }
\keyword{ pca }% __ONLY ONE__ keyword per line
