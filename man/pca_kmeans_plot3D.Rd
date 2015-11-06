\name{pca_kmeans_plot3D}
\alias{pca_kmeans_plot3D}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  3D PCA k-means plot (interactive)
}
\description{
  Groups the points with the clusters given by k-means in a interactive 3D PCA scores plot.
}
\usage{
pca_kmeans_plot3D(dataset, pca.result, num.clusters = 3, 
pcas = c(1, 2, 3), kmeans.result = NULL, labels = FALSE, 
size = 1, ...)
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
  \item{size}{
    parameter of plot3d from rgl package.
}
  \item{\dots}{
	additional parameters of plot3d function from rgl package.
}
}

\examples{
  ## Example of a 3D PCA k-means plot
  data(cachexia)
  pca.result = pca_analysis_dataset(cachexia)
  pca_kmeans_plot3D(cachexia, pca.result, num.clusters = 3, 
    		    pcas = c(1,2,3))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ kmeans }
\keyword{ pca }% __ONLY ONE__ keyword per line
