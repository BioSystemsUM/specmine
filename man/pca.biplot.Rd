\name{pca.biplot}
\alias{pca.biplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  PCA biplot
}
\description{
  Shows a PCA biplot.
}
\usage{
pca.biplot(pca.result, cex = 0.8, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pca.result}{
            prcomp object with the PCA results.

}
  \item{cex}{
    value that indicates the magnification relative to the default.
}
  \item{\dots}{
    additional parameters passed to biplot function.
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
  ## Example of a PCA biplot
  pca.biplot(pca.result, cex = 0.8)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~pca }
\keyword{ ~biplot }% __ONLY ONE__ keyword per line
