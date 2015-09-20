\name{absorbance.to.transmittance}
\alias{absorbance.to.transmittance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
   Convert absorbance to transmittance
}
\description{
   Converts absorbance values to transmittance values.
}
\usage{
absorbance.to.transmittance(dataset)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
}

\value{
   Returns the dataset with the data points converted to transmittance values.
}

\examples{
\dontrun{
  ## Example of converting transmittance values to absorbance values
  dataset = absorbance.to.transmittance(dataset)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~absorbance }
\keyword{ ~transmittance }% __ONLY ONE__ keyword per line
