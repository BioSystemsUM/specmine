\name{dendrogram.plot}
\alias{dendrogram.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Plot dendrogram
}
\description{
  Plot dendrogram of hierarchical clustering results.
}
\usage{
dendrogram.plot(dataset, hc.result, column.metadata = 1, labels = NULL, ...)
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
  ## Example of a dendrogram
  dendrogram.plot(metabolomics.experiment, hc.result)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~clustering }
\keyword{ ~dendrogram }
\keyword{ ~hclust }% __ONLY ONE__ keyword per line
