\name{dendrogram.plot.col}
\alias{dendrogram.plot.col}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Plot dendrogram
}
\description{
  Plot dendrogram of hierarchical clustering results with different colors
}
\usage{
dendrogram.plot.col(dataset, hc.result, classes.col, title = "", lab.cex = 1, leg.pos = "topright",...)
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
    string or index indicating what metadata to use to name the leafs.
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
  \item{\dots}{
    other parameters for plotting.
}
}

\examples{
\dontrun{
  ## Example of colored dendrogram
  dendrogram.plot.col(metabolomics.dataset, hc.result, "class", title = "Example")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~clustering }
\keyword{ ~dendrogram }
\keyword{ ~hclust }
