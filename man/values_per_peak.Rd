\name{values_per_peak}
\alias{values_per_peak}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Values per peak
}
\description{
  Gets the number of values on each peak.
}
\usage{
values_per_peak(samples.df)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{samples.df}{
    data frame with the samples' peaks.
}
}

\value{
  Returns a vector with the number of values for each peak.
}

\examples{
  ## Example of getting the number of values for each peak
  library(specmine.datasets)
  data(propolis)
  num.peaks = values_per_peak(propolis$data)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ peak }
\keyword{ sample }% __ONLY ONE__ keyword per line
