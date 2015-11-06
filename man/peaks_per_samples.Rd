\name{peaks_per_samples}
\alias{peaks_per_samples}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Peaks per samples
}
\description{
  Calculates the number of peaks on each sample.
}
\usage{
peaks_per_samples(sample.list)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sample.list}{
    list of data frames with the samples' peaks.
}
}

\value{
  Returns a numeric vector with the number of peaks on each sample.
}

\examples{
    ## Example of counting the peaks in each sample
    data(propolisSampleList)
    num.peaks.samples = peaks_per_samples(propolisSampleList)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ peak }
\keyword{ sample }% __ONLY ONE__ keyword per line
