\name{hierarchical_clustering}
\alias{hierarchical_clustering}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Perform hierarchical clustering analysis
}
\description{
  Perform hierarchical clustering analysis on the dataset.
}
\usage{
hierarchical_clustering(dataset, distance = "euclidean", 
clustMethod = "complete", hc.type = "samples")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{distance}{
    the distance measure to be used to compute the distances between the rows of a data matrix. Possible types are "euclidean", "manhattan", "pearson" or "spearman".
}
  \item{clustMethod}{
    the agglomeration method to be used. Possible values are "ward", "single", "complete", "average", "mcquitty", "median" or "centroid".
}
  \item{hc.type}{
    a string indicating if hierarchical cluster analysis will be performed on samples ("samples") or on variables ("variables")
}
}

\value{
  An object of class hclust with the clustering results.
}

\examples{
  ## Example of hierarchical clustering
  data(cachexia)
  hc.result = hierarchical_clustering(cachexia, 
	      distance = "euclidean", clustMethod = "complete", 
	      hc.type = "samples")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ clustering }
\keyword{ hclust }% __ONLY ONE__ keyword per line
