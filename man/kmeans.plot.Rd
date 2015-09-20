\name{kmeans.plot}
\alias{kmeans.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Plot kmeans clusters
}
\description{
  Plot for each formed cluster, in grey the values of all samples of that cluster and in blue the median of that samples.
}
\usage{
kmeans.plot(dataset, kmeans.result)
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
\dontrun{
  ## Example of kmeans plot
  kmeans.plot(metabolomics.dataset, kmeans.result)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~clustering }
\keyword{ ~kmeans }% __ONLY ONE__ keyword per line
\keyword{ ~plotting }