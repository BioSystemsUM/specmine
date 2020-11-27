\name{aggregate_samples}
\alias{aggregate_samples}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Aggregate samples
}
\description{
 Aggregate samples according to an aggregate function like mean, median, etc.
}
\usage{
aggregate_samples(dataset, indexes, aggreg.fn = "mean", 
meta.to.remove = c())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{indexes}{
    index vector with the samples that are going to be aggregated (e.g. c(1,1,2,2), this index vector will aggregate the first two samples and the last two samples).
}
  \item{aggreg.fn}{
    aggregation function (e.g. "mean", "median", etc).
}
  \item{meta.to.remove}{
    metadata's variables to be removed.
}
}

\value{
 Returns the dataset with the samples aggregated.
}


\examples{
  ## Example of aggregating samples
  library(specmine.datasets)
  data(propolis)
  propolis_proc = missingvalues_imputation(propolis)
  dataset = aggregate_samples(propolis_proc, as.integer(propolis$metadata$seasons), "mean") 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ aggregation }
\keyword{ sample }% __ONLY ONE__ keyword per line
