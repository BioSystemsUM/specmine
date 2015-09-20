\name{is.spectra}
\alias{is.spectra}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Check type of data
}
\description{
  Check if the dataset is from spectral data where x.values are numeric.
}
\usage{
is.spectra(dataset)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
}

\value{
  Returns a boolean indicating if the dataset is from spectral data or not.
}

\examples{
\dontrun{
  ## Example of checking if the dataset is from spectral data
  is.spectra(dataset)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~dataset }
\keyword{ ~spectral }% __ONLY ONE__ keyword per line
