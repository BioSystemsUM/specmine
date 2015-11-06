\name{dendrogram_plot}
\alias{dendrogram_plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Plot dendrogram
}
\description{
  Plot dendrogram of hierarchical clustering results.
}
\usage{
dendrogram_plot(dataset, hc.result, column.metadata = 1, 
labels = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{hc.result}{
    object of class hclust with the clustering results.
}
  \item{column.metadata}{
    string or index indicating what metadata to use to name the leafs.
}
  \item{labels}{
    vector with the leaf names (optional).
}
  \item{\dots}{
    other parameters for plotting.
}
}

\examples{
  ## Example of a dendrogram
  data(cachexia)
  hc.result = hierarchical_clustering(cachexia)
  dendrogram_plot(cachexia, hc.result)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ clustering }
\keyword{ dendrogram }
\keyword{ hclust }% __ONLY ONE__ keyword per line
