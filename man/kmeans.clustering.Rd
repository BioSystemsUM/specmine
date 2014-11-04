\name{kmeans.clustering}
\alias{kmeans.clustering}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Perform k-means clustering analysis
}
\description{
Perform k-means clustering analysis on the dataset.
}
\usage{
kmeans.clustering(dataset, num.clusters, type = "samples")
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
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  An object of class kmeans with the clustering results.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
\dontrun{
  ## Example of kmeans clustering
  kmeans.result = kmeans.clustering(metabolomics.dataset, num.clusters = 4, type = "samples")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~clustering }
\keyword{ ~kmeans }% __ONLY ONE__ keyword per line
