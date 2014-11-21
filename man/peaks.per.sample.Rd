\name{peaks.per.sample}
\alias{peaks.per.sample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Peaks per sample
}
\description{
  Counts number of peaks in a sample (given its index).
}
\usage{
peaks.per.sample(sample.list, sample.index)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sample.list}{
    list of data frames with the samples' peaks.
}
  \item{sample.index}{
    sample index.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  Return a integer value with the number of peaks in the sample.
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
    ## Example of counting the peaks in a sample
    num.peaks.sample = peaks.per.sample(sample.list, 4)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~peak }
\keyword{ ~sample }% __ONLY ONE__ keyword per line
