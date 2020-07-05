\name{kmeans_result_df}
\alias{kmeans_result_df}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Show cluster's members  
}
\description{
  Show for each cluster from kmeans analysis the sample names belonging to them.
}
\usage{
kmeans_result_df(kmeans.result)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{kmeans.result}{
    object of class kmeans with the clustering results.
}
}

\value{
Data frame with the clusters and the samples' names that belong to each one.
}

\examples{
  ## Example of showing kmeans cluster's members
  kmeans.result = kmeans_clustering(cachexia, 
		  num.clusters = 4, type = "samples")
  kmeans_result_df(kmeans.result)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ clustering }
\keyword{ kmeans }% __ONLY ONE__ keyword per line
