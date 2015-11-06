\name{num_samples}
\alias{num_samples}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Get number of samples
}
\description{
  Get the number of samples from a dataset.
}
\usage{
num_samples(dataset)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
}

\value{
  Returns an integer with the number of samples in the dataset.
}

\examples{
  ## Example of getting the number of samples
  data(cachexia)
  number.of.samples = num_samples(cachexia)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ samples }
\keyword{ dataset }% __ONLY ONE__ keyword per line
