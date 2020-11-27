\name{dendrogram_plot_col}
\alias{dendrogram_plot_col}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Plot dendrogram
}
\description{
  Plot dendrogram of hierarchical clustering results with different colors
}
\usage{
dendrogram_plot_col(dataset, hc.result, classes.col, colors = NULL, title = "", 
lab.cex = 1, leg.pos = "topright", label_samples=NULL,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{hc.result}{
    object of class hclust with the clustering results.
}
  \item{classes.col}{
    string or index indicating what metadata to use to color the leafs.
}
  \item{colors}{
    vector with the corresponding colors of the metadata classes.
}
  \item{title}{
    title of dendrogram.    
}
  \item{lab.cex}{
    the magnification to be used for x and y labels relative to the current setting of cex.
}
  \item{leg.pos}{
    position of the legend.
}
  \item{label_samples}{
    string or index indicating what metadata to use to name the leafs. If not provided the name of the leafs will remain the samples names.
}
  \item{\dots}{
    other parameters for plotting.
}
}

\examples{
  ## Example of colored dendrogram
  library(specmine.datasets)
  data(cachexia)
  hc.result = hierarchical_clustering(cachexia)
  dendrogram_plot_col(cachexia, hc.result, "Muscle.loss", 
  title = "Example")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ clustering }
\keyword{ dendrogram }
\keyword{ hclust }
