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
clustering(dataset, method = "hc", distance = "euclidean", type = "samples", num.clusters = 5)
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
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  An object of class kmeans or hclust with the clustering results.
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
  ## Example of kmeans and hierarchical clustering
  hc.result = clustering(metabolomics.dataset, method = "hc", distance = "euclidean")
  kmeans.result = clustering(metabolomics.dataset, method = "kmeans", num.clusters = 4)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~clustering }
\keyword{ ~kmeans }% __ONLY ONE__ keyword per line
\keyword{ ~hclust }