\name{savitzky.golay}
\alias{savitzky.golay}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Savitzky-golay transformation
}
\description{
  Smoothing and derivative of the data using Savitzky-Golay.
}
\usage{
savitzky.golay(dataset, p.order, window, deriv = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
  list representing the dataset from a metabolomics experiment.
}
  \item{p.order}{
  integer value representing the polynomial order.
}
  \item{window}{
  odd number indicating the window size.
}
  \item{deriv}{
  integer value indicating the differentiation order.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns the dataset with the spectra smoothed using Savitzky-Golay.
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
  ## Example of smoothing the spectra from dataset with Savitzky-Golay
  dataset.smoothed = savitzky.golay(dataset, p.order = 3, window = 11, deriv = 0)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~smoothing }
\keyword{ ~savitzky-golay }% __ONLY ONE__ keyword per line
