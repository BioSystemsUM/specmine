\name{savitzky_golay}
\alias{savitzky_golay}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Savitzky-golay transformation
}
\description{
  Smoothing and derivative of the data using Savitzky-Golay.
}
\usage{
savitzky_golay(dataset, p.order, window, deriv = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
  list representing the dataset from a metabolomics experiment.
}
  \item{p.order}{
  integer value representing the polynomial order.
}
  \item{window}{
  odd number indicating the window size.
}
  \item{deriv}{
  integer value indicating the differentiation order.
}
}

\value{
Returns the dataset with the spectra smoothed using Savitzky-Golay.
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ smoothing }
\keyword{ savitzky-golay }% __ONLY ONE__ keyword per line
