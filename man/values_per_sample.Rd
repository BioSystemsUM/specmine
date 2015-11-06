\name{values_per_sample}
\alias{values_per_sample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Values per peak
}
\description{
  Gets the number of values on each sample.
}
\usage{
values_per_sample(samples.df)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{samples.df}{
    data frame with the samples' peaks.
}
}

\value{
  Returns a vector with the number of values for each sample.
}

\examples{
  ## Example of getting the number of values for each sample
  data(propolis)
  num.values.samples = values_per_sample(propolis$data)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ values }
\keyword{ sample }% __ONLY ONE__ keyword per line
