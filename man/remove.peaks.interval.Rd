\name{remove.peaks.interval}
\alias{remove.peaks.interval}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Remove interval of peaks
}
\description{
  Removes peaks from a given interval.
}
\usage{
remove.peaks.interval(sample.df, peak.val.min, peak.val.max)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sample.df}{
    data frame with the samples' peaks.
}
  \item{peak.val.min}{
    minimum peak value.
}
  \item{peak.val.max}{
    maximum peak value.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  Returns a data frame with the filtered samples' peaks. 
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
  ## Example of removing a interval of peaks
  samples.df = remove.peaks.interval(sample.df, 2, 8.05)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~peak }
\keyword{ ~sample }% __ONLY ONE__ keyword per line
