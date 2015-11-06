\name{absorbance_to_transmittance}
\alias{absorbance_to_transmittance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
   Convert absorbance to transmittance
}
\description{
   Converts absorbance values to transmittance values.
}
\usage{
absorbance_to_transmittance(dataset)
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
  ## Example of converting transmittance values to absorbance values
  data(cassavaPPD)
  cassavaPPD = transmittance_to_absorbance(cassavaPPD)
  cassavaPPD = absorbance_to_transmittance(cassavaPPD)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ absorbance }
\keyword{ transmittance }% __ONLY ONE__ keyword per line
