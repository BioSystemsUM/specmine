\name{clustering}
\alias{clustering}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Perform cluster analysis
}
\description{
Perform cluster analysis on the dataset.
}
\usage{
clustering(dataset, method = "hc", distance = "euclidean", 
type = "samples", num.clusters = 5, clustMethod = "complete")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{method}{
    a string describing the method of clustering. Possible types are:
    \itemize{
      \item "hc" perform hierarchical clustering.
      \item "kmeans" perform kmeans clustering.
    }
}
  \item{distance}{
    the distance measure to be used to compute the distances between the rows of a data matrix. Possible types are "euclidean", "manhattan", "pearson" or "spearman". Only for hierarchical clustering.
}
  \item{type}{
    a string indicating if cluster analysis will be performed on samples ("samples") or on variables ("variables"). 
}
  \item{num.clusters}{
    the number of clusters in k-means cluster analysis.
}
  \item{clustMethod}{
    Cluster method for hierarchical clustering.
}
}

\value{
  An object of class kmeans or hclust with the clustering results.
}

\examples{
  ## Example of kmeans and hierarchical clustering
  library(specmine.datasets)
  data(cachexia)
  hc.result = clustering(cachexia, method = "hc", 
	      distance = "euclidean")
  kmeans.result = clustering(cachexia, method = "kmeans", 
		  num.clusters = 4)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ clustering }
\keyword{ kmeans }% __ONLY ONE__ keyword per line
\keyword{ hclust }
