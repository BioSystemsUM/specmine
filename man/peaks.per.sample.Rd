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

\value{
  Return a integer value with the number of peaks in the sample.
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
