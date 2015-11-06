\name{kmeans_clustering}
\alias{kmeans_clustering}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Perform k-means clustering analysis
}
\description{
Perform k-means clustering analysis on the dataset.
}
\usage{
kmeans_clustering(dataset, num.clusters, type = "samples")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{num.clusters}{
    the number of clusters.
}
  \item{type}{   
    a string indicating if k-means will be performed on samples ("samples") or on variables ("variables")
}
}

\value{
  An object of class kmeans with the clustering results.
}

\examples{
  ## Example of kmeans clustering
  data(cachexia)
  kmeans.result = kmeans_clustering(cachexia, 
		  num.clusters = 4, type = "samples")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ clustering }
\keyword{ kmeans }% __ONLY ONE__ keyword per line
