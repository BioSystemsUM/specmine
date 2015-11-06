\name{get_peak_values}
\alias{get_peak_values}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Get peak values
}
\description{
  Gets the peak values from a data frame of samples' peaks.
}
\usage{
get_peak_values(samples.df, peak.val)
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

\value{
  Returns a vector with the peak values.
}

\examples{
  ## Example of getting the peak values
  data(propolis)
  peak.values = get_peak_values(propolis$data, 2.11)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ peak }
\keyword{ sample }% __ONLY ONE__ keyword per line
