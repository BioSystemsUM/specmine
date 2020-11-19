\name{smoothing_interpolation}
\alias{smoothing_interpolation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Smoothing interpolation
}
\description{
  Performs smoothing interpolation according to a specific method.
}
\usage{
smoothing_interpolation(dataset, method = "bin", 
reducing.factor = 2, x.axis = NULL, p.order = 3, 
window = 11, deriv = 0, na.rm = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
  list representing the dataset from a metabolomics experiment.
}
  \item{method}{
  string specifying the smoothing method. The three possible methods are: "bin", "loess" and "savitzky.golay".
}
  \item{reducing.factor}{
  numeric value indicating the reducing factor for bin method.
}
  \item{x.axis}{
  numeric vector representing the new x-axis for loess method.
}
  \item{p.order}{
  numeric value representing the polynomial order for savitzky-golay method. 
  }
  \item{window}{
  numeric value indicating the size of the window for savitzky-golay method. (must be an odd number)
  }
  \item{deriv}{
  numeric value indicating the differentiation order for savitzky-golay method. 
  }
  \item{na.rm}{
  boolean value indicating if \code{NA}s should be removed. Defaults to TRUE.
  }
}

\value{
  Returns the dataset with the spectra smoothed.
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ smoothing }
\keyword{ bin }% __ONLY ONE__ keyword per line
\keyword{ loess }
\keyword{ savitzky-golay }
