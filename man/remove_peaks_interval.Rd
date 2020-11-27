\name{remove_peaks_interval}
\alias{remove_peaks_interval}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Remove interval of peaks
}
\description{
  Removes peaks from a given interval.
}
\usage{
remove_peaks_interval(sample.df, peak.val.min, peak.val.max)
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

\value{
  Returns a data frame with the filtered samples' peaks. 
}

\examples{
  ## Example of removing a interval of peaks
  library(specmine.datasets)
  data(propolisSampleList)
  samples.df = remove_peaks_interval(propolisSampleList[[1]], 2, 8.05)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ peak }
\keyword{ sample }% __ONLY ONE__ keyword per line
