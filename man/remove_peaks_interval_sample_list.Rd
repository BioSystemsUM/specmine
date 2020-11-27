\name{remove_peaks_interval_sample_list}
\alias{remove_peaks_interval_sample_list}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Remove interval of peaks (sample list)
}
\description{
  Removes peaks on a sample list given a peak interval.
}
\usage{
remove_peaks_interval_sample_list(sample.list, peak.val.min, 
peak.val.max)
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

\value{
  Returns the sample list with the filtered peaks.
}

\examples{
  ## Example of removing a interval of peaks in a sample list
  library(specmine.datasets)
  data(propolisSampleList)
  samples.list = remove_peaks_interval_sample_list(propolisSampleList, 2, 8.05)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ peak }
\keyword{ sample }% __ONLY ONE__ keyword per line
