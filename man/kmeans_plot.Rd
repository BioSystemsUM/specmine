\name{kmeans_plot}
\alias{kmeans_plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Plot kmeans clusters
}
\description{
  Plot for each formed cluster, in grey the values of all samples of that cluster and in blue the median of that samples.
}
\usage{
kmeans_plot(dataset, kmeans.result)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{kmeans.result}{
    object of class kmeans with the clustering results.
}
}


\examples{
  ## Example of kmeans plot - dataset filtered for performance purposes
  library(specmine.datasets)
  data(cachexia)
  kmeans.result = kmeans_clustering(cachexia, 
		  num.clusters = 4, type = "samples")
  kmeans_plot(cachexia, kmeans.result)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ clustering }
\keyword{ kmeans }% __ONLY ONE__ keyword per line
\keyword{ plotting }
