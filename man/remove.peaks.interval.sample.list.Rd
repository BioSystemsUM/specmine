\name{remove.peaks.interval.sample.list}
\alias{remove.peaks.interval.sample.list}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Remove interval of peaks (sample list)
}
\description{
  Removes peaks on a sample list given a peak interval.
}
\usage{
remove.peaks.interval.sample.list(sample.list, peak.val.min, peak.val.max)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sample.list}{
    list of data frames with the samples' peaks.
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
  Returns the sample list with the filtered peaks.
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
  ## Example of removing a interval of peaks in a sample list
  samples.df = remove.peaks.interval.sample.list(sample.list, 2, 8.05)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~peak }
\keyword{ ~sample }% __ONLY ONE__ keyword per line
