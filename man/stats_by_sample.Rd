\name{stats_by_sample}
\alias{stats_by_sample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Statistics of samples
}
\description{
  Get a summary of statistics of the samples.
}
\usage{
stats_by_sample(dataset, samples = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{samples}{
    if defined restricts the application to a given set of samples.
}
}

\value{
  Returns a vector with the a summary of statistics of the samples.
}

\examples{
   \dontrun{
     ## Example of getting stats of samples
     samples.stats.result = stats_by_samples(dataset)
   }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ statistic }
\keyword{ sample }% __ONLY ONE__ keyword per line
