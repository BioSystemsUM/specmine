\name{pca.kmeans.plot2D}
\alias{pca.kmeans.plot2D}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  2D PCA k-means plot
}
\description{
    Groups the points with the clusters given by k-means in a 2D PCA scores plot.
}
\usage{
pca.kmeans.plot2D(dataset, pca.result, num.clusters = 3, pcas = c(1, 2), kmeans.result = NULL, labels = FALSE, ellipses = FALSE)
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
  \item{pcas}{
    vector with the principal components to be plotted.
}
  \item{kmeans.result}{
    result from k-means. If null k-means is performed in the function.
}
  \item{labels}{
    boolean value indicating if the samples' labels will be shown.
}
  \item{ellipses}{
    boolean value that indicates if an ellipse will be drawn on each group of the metadata's variable.
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
  ## Example of a 2D PCA k-means plot
  pca.kmeans.plot2D(dataset, pca.result, num.clusters = 3, pcas = c(1,2))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~pca }
\keyword{ ~kmeans }% __ONLY ONE__ keyword per line
