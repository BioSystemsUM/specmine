\name{transmittance_to_absorbance}
\alias{transmittance_to_absorbance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
   Convert transmittance to absorbance
}
\description{
   Converts transmittance values to absorbance values.
}
\usage{
transmittance_to_absorbance(dataset, percent = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{percent}{
    boolean value to indicate if the absorbance values are in a percentage or not.
}
}

\value{
   Returns the dataset with the data points converted to absorbance values.
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ transmittance }
\keyword{ absorbance }% __ONLY ONE__ keyword per line
