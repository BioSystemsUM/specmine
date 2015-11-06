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
  data(cassavaPPD)
  index.vector = c(1,1,2,2,2,2,2,1,1,3,3,3,3,3,1,4,4,4,4,4,5,5,6,6,6,6,6,5,5,
		  7,7,7,7,7,8,8,8,8,8,9,9,10,10,10,10,10,9,9,11,11,11,11,11,
		  9,12,12,12,12,12,13,13,14,14,14,14,14,13,13,15,15,15,15,15,
		  13,16,16,16,16,16)
  dataset = aggregate_samples(cassavaPPD, index.vector, "mean") 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ aggregation }
\keyword{ sample }% __ONLY ONE__ keyword per line
