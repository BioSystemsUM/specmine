\name{plot.spectra}
\alias{plot.spectra}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot spectra
}
\description{
Plot spectra from dataset.
}
\usage{
plot.spectra(dataset, column.class, func = NULL, samples = NULL, variable.bounds = NULL, xlab = NULL, ylab = NULL, lty = 1, legend.place = "topright", cex = 0.0, reverse.x = F, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{column.class}{
    string indicating the metadata's variable.
}
  \item{func}{
    function to compute the summary statistics to apply to the data.
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
    y-axis label.
}
  \item{lty}{
    parameter of matplot.
}
  \item{legend.place}{
    string indicating the place that the legend's box will be placed.
}
  \item{cex}{
    numeric value that indicates the amount by which the legend is magnified relative to the default.
}
  \item{reverse.x}{
    boolean value indicating if the x-axis will be shown reversed or not.  
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
  ## Example of plotting spectra (simple)
  plot.spectra(dataset, "metadata.var1", func = NULL, samples = c("sample1", "sample2"), variable.bounds = c(200,2000), legend.place = "topright")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~spectra }
\keyword{ ~plot }% __ONLY ONE__ keyword per line
