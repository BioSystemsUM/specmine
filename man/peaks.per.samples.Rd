\name{peaks.per.samples}
\alias{peaks.per.samples}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Peaks per samples
}
\description{
  Calculates the number of peaks on each sample.
}
\usage{
peaks.per.samples(sample.list)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sample.list}{
    list of data frames with the samples' peaks.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  Returns a numeric vector with the number of peaks on each sample.
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
    ## Example of counting the peaks in each sample
    num.peaks.samples = peaks.per.samples(sample.list)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~peak }
\keyword{ ~sample }% __ONLY ONE__ keyword per line
