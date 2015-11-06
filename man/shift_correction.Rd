\name{shift_correction}
\alias{shift_correction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Shift correction
}
\description{
  Shifts the spectra according to a specific method.
}
\usage{
shift_correction(dataset, method = "constant", shift.val = 0, 
interp.function = "linear", ref.limits = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{method}{
    string that indicates the shifting method. It can be:
    \itemize{
      \item "constant" uses a constant shift that is added to the x-values
      \item "interpolation" uses interpolation according to "interp.function"
    }
}
  \item{shift.val}{
     value of the shift (for constant and interpolation methods); can be a single value for all spectra, can be the string "auto", the shifts are automatically determined or a vector with the size of the number of samples with the shifts for each spectra.
}
  \item{interp.function}{
  string that represents the interpolation function, can be "linear" or "spline".
}
  \item{ref.limits}{
  vector with 2 elements that represents the reference limits to calculate the shifts.
}
}

\value{
Returns the dataset with the spectras shifted.
}

\examples{
  ## Example of shift correction with linear interpolation and 
  ## shifts calculated based on a reference band of the spectra
  data(cassavaPPD)
  dataset.corrected = shift_correction(cassavaPPD, method = "interpolation", 
		      shift.val = "auto", interp.function = "linear", 
		      ref.limits = c(800,850))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ shift }
\keyword{ correction }% __ONLY ONE__ keyword per line
