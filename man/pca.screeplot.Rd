\name{pca.screeplot}
\alias{pca.screeplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  PCA scree plot
}
\description{
  PCA scree plot with the proportion and cumulative variance of the PCs.
}
\usage{
pca.screeplot(pca.result, num.pcs = NULL, cex.leg = 0.8, leg.pos = "right", lab.text = c("individual percent", "cumulative percent"), fill.col = c("blue", "red"), ylab = "Percentage", xlab = "Principal components", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pca.result}{
    prcomp object with the PCA results.
}
  \item{num.pcs}{
    number of principal components.
}
  \item{cex.leg}{
    cex value of legend.
}
  \item{leg.pos}{
    legend position.
}
  \item{lab.text}{
    legend's labels.
}
  \item{fill.col}{
    color of the legend's boxes.
}
  \item{ylab}{
    y-axis label.
}
  \item{xlab}{
    x-axis label
}
  \item{\dots}{
    additional parameters to matplot.
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
  ## Example of a scree plot
  pca.screeplot(pca.result)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~screeplot }
\keyword{ ~pca }% __ONLY ONE__ keyword per line
