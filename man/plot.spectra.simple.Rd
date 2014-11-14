\name{plot.spectra.simple}
\alias{plot.spectra.simple}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Plot spectra (simple)
}
\description{
  Plot spectra from dataset (simple version).
}
\usage{
plot.spectra.simple(dataset, samples = NULL, variable.bounds = NULL, xlab = NULL, ylab = NULL, lty = 1, lwd = 1, col = 1, reverse.x = F, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{samples}{
    vector with samples' names, if NULL all the samples will be considered.
}
  \item{variable.bounds}{
    numeric vector with two elements indicating the interval of x-values to plot.
}
  \item{xlab}{
    x-axis label.
}
  \item{ylab}{
    y-axis label
}
  \item{lty}{
    parameter of matplot.
}
  \item{lwd}{
    parameter of matplot.
}
  \item{col}{
    parameter of matplot.
}
  \item(reverse.x){
    boolean value indicating if the x-axis will be shown reversed or not.  
}
  \item{\dots}{
    additional parameters to pass to matplot.
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
  ## Example of plotting spectra (simple)
  plot.spectra.simple(dataset, samples = c("sample1", "sample2"), variable.bounds = c(200,2000))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~plot }
\keyword{ ~spectra }% __ONLY ONE__ keyword per line
