\name{pca.biplot3D}
\alias{pca.biplot3D}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  3D PCA biplot (interactive)
}
\description{
  Shows a interactive 3D PCA biplot.
}
\usage{
pca.biplot3D(dataset, pca.result, column.class, pcas = c(1, 2, 3))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
        list representing the dataset from a metabolomics experiment.
}
  \item{pca.result}{
    prcomp object with the PCA results.
}
  \item{column.class}{
    metadata's variable.
}
  \item{pcas}{
    the three principal components.
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
  ## Example of a 3D PCA biplot
  pca.biplot3D(dataset, pca.result, "metadata.var2", pcas = c(1,2,3))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~pca }
\keyword{ ~biplot }% __ONLY ONE__ keyword per line
