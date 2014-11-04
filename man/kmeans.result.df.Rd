\name{kmeans.result.df}
\alias{kmeans.result.df}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Show cluster's members  
}
\description{
  Show for each cluster from kmeans analysis the sample names belonging to them.
}
\usage{
kmeans.result.df(kmeans.result, centers)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{kmeans.result}{
    object of class kmeans with the clustering results.
}
  \item{centers}{
    number of clusters
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Data frame with the clusters and the samples' names that belong to each one.
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
  ## Example of showing kmeans cluster's members
  kmeans.result.df = function(kmeans.result, 4)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~clustering }
\keyword{ ~kmeans }% __ONLY ONE__ keyword per line
