\name{pca.pairs.kmeans.plot}
\alias{pca.pairs.kmeans.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  PCA k-means pairs plot
}
\description{
  Groups the points with the clusters from k-means in a PCA pairs plot.
}
\usage{
pca.pairs.kmeans.plot(dataset, pca.result, num.clusters = 3, kmeans.result = NULL, pcas = c(1, 2, 3, 4, 5))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{pca.result}{
    prcomp object with the PCA results.
}
  \item{num.clusters}{
    number of clusters of k-means.
}
  \item{kmeans.result}{
    result from k-means. If null k-means is performed in the function.
}
  \item{pcas}{
    vector with the principal components to be plotted.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
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
  ## Example of a PCA k-means pairs plot
  pca.pairs.kmeans.plot(dataset, pca.result, num.clusters = 3, kmeans.result = kmeans.res, pcas = c(1,2,3,4,5))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kmeans }
\keyword{ ~pairs }% __ONLY ONE__ keyword per line
