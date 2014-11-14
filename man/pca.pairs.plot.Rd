\name{pca.pairs.plot}
\alias{pca.pairs.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  PCA pairs plot
}
\description{
  Shows a PCA pairs plot.
}
\usage{
pca.pairs.plot(dataset, column.class, pca.result, pcas = c(1, 2, 3, 4, 5))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
        list representing the dataset from a metabolomics experiment.
}
  \item{column.class}{
    metadata's variable.
}
  \item{pca.result}{
    prcomp object with the PCA results.
}
  \item{pcas}{
    the principal components to be shown.
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
  ## Example of a PCA pairs plot
  pca.pairs.plot(dataset, "metadata.var2", pca.result, pcas = c(1,2,3))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~pairs }
\keyword{ ~pca }% __ONLY ONE__ keyword per line
