\name{smoothing.spcbin.hyperspec}
\alias{smoothing.spcbin.hyperspec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Wavelength binning
}
\description{
Smoothing interpolation with wavelength binning.
}
\usage{
smoothing.spcbin.hyperspec(dataset, reducing.factor = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
  list representing the dataset from a metabolomics experiment.
}
  \item{reducing.factor}{
  numerical value indicating the reducing factor.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  Returns the dataset with the spectra smoothed with the binning method.
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
  ## Example of smoothing the spectra from dataset with binning
  dataset.smoothed = smoothing.spcbin.hyperspec(dataset, reducing.factor = 4)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~binning }
\keyword{ ~smoothing }% __ONLY ONE__ keyword per line
