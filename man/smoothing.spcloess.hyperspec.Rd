\name{smoothing.spcloess.hyperspec}
\alias{smoothing.spcloess.hyperspec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Loess smoothing
}
\description{
  Smoothing interpolation on a new wavelength axis using loess.
}
\usage{
smoothing.spcloess.hyperspec(dataset, x.axis = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
  list representing the dataset from a metabolomics experiment.
}
  \item{x.axis}{
  numeric vector that represents the new wavelength axis to interpolate on.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns the dataset with the spectra smoothed using loess.
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
  ## Example of smoothing the spectra from dataset with loess
  new.axis = seq(200,800,by=4)
  dataset.smoothed = smoothing.spcloess.hyperspec(dataset, x.axis = new.axis)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~smoothing }
\keyword{ ~loess }% __ONLY ONE__ keyword per line
