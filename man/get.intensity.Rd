\name{get.intensity}
\alias{get.intensity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Get intensity
}
\description{
  Get the value of an internsity given the sample (data frame) and the frequency.
}
\usage{
get.intensity(sample.df, freq, tolerance = 0.001)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sample.df}{
    data frame with the samples' peaks.
}
  \item{freq}{
    frequence.
  }
  \item{tolerance}{
    tolerance value of the frequency.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
    Returns the intensity (or a NA value in case that the frequency doesn't exist).
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
    ## Example of getting the intensity
    frequencies = get.intensity(sample, 2.01, tolerance = 0.001)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~frequencies }
\keyword{ ~intensity }% __ONLY ONE__ keyword per line
