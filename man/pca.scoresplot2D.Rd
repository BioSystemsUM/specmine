\name{pca.scoresplot2D}
\alias{pca.scoresplot2D}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  2D PCA scores plot
}
\description{
  Shows a 2D PCA scores plot of two principal componets.
}
\usage{
pca.scoresplot2D(dataset, pca.result, column.class, pcas = c(1, 2), labels = FALSE, ellipses = FALSE, pallette = 2)
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
    vector of two elements with the PCs that will be plotted.
}
  \item{labels}{
    boolean value indicating if the sample's labels will be displayed.
}
  \item{ellipses}{
    boolean value that indicates if an ellipse will be drawn on each group of the metadata's variable.
}
  \item{pallette}{
    parameter of scale_colour_brewer from ggplot2.
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
  ## Example of a 2D PCA scores plot
  pca.scoresplot2D(dataset, pca.result, "metadata.var1", pcas = c(1,2), ellipses = T)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~scoresplot }
\keyword{ ~pca }% __ONLY ONE__ keyword per line
