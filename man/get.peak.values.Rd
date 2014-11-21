\name{get.peak.values}
\alias{get.peak.values}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Get peak values
}
\description{
  Gets the peak values from a data frame of samples' peaks.
}
\usage{
get.peak.values(samples.df, peak.val)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{samples.df}{
    data frame with the samples' peaks.
}
  \item{peak.val}{
    peak name.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  Returns a vector with the peak values.
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
  ## Example of getting the peak values
  peak.values = get.peak.values(samples.df, 2.10)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~peak }
\keyword{ ~sample }% __ONLY ONE__ keyword per line
