\name{pca_kmeans_plot2D}
\alias{pca_kmeans_plot2D}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  2D PCA k-means plot
}
\description{
    Groups the points with the clusters given by k-means in a 2D PCA scores plot.
}
\usage{
pca_kmeans_plot2D(dataset, pca.result, num.clusters = 3, 
pcas = c(1, 2), kmeans.result = NULL, labels = FALSE, bw=F,
ellipses = FALSE, leg.pos = "right", xlim = NULL, ylim = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{pca.result}{
    prcomp object with the PCA results.
}
  \item{num.clusters}{
    number of clusters of k-means.
}
  \item{pcas}{
    vector with the principal components to be plotted.
}
  \item{kmeans.result}{
    result from k-means. If null k-means is performed in the function.
}
  \item{labels}{
    boolean value indicating if the samples' labels will be shown.
}
  \item{ellipses}{
    boolean value that indicates if an ellipse will be drawn on each group of the metadata's variable. Ellipses will not be drawn if bw=TRUE.
}
  \item{bw}{
    if TRUE, it will be displayed a black and white plot. It defaults to FALSE.
}
  \item{leg.pos}{
    legend position.
}
  \item{xlim}{
    vector with two positions with the x-axis limits.
}
  \item{ylim}{
    vector with two positions with the y-axis limits.
}
}

\examples{
\dontrun{
  ## Example of a 2D PCA k-means plot
  data(cachexia)
  pca.result = pca_analysis_dataset(cachexia)
  pca_kmeans_plot2D(cachexia, pca.result, num.clusters = 3, pcas = c(1,2))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ pca }
\keyword{ kmeans }% __ONLY ONE__ keyword per line
