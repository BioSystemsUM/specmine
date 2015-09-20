\name{num.samples}
\alias{num.samples}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Get number of samples
}
\description{
  Get the number of samples from a dataset.
}
\usage{
num.samples(dataset)
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
\dontrun{
  ## Example of getting the number of samples
  number.of.samples = num.samples(dataset)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~samples }
\keyword{ ~dataset }% __ONLY ONE__ keyword per line
