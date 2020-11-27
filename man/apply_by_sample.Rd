\name{apply_by_sample}
\alias{apply_by_sample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Apply function to samples
}
\description{
  Applies a function to the values of each sample
}
\usage{
apply_by_sample(dataset, fn.to.apply, samples = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{fn.to.apply}{
    function to apply (e.g. mean, max, min).
}
  \item{samples}{
    if defined restricts the application to a given set of samples.
}
  \item{\dots}{
    additional parameters to apply function.
}
}


\value{
   Returns a vector with the samples and the value of the applied function.
}


\examples{
     ## Example of applying a function to variables
	 library(specmine.datasets)
     data(cachexia)
     apply.samples.result = apply_by_sample(cachexia, mean)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ apply }
\keyword{ sample }% __ONLY ONE__ keyword per line
