\name{baseline_correction}
\alias{baseline_correction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Baseline correction
}
\description{
  Performs baseline correction on the dataset.
}
\usage{
baseline_correction(dataset, method = "modpolyfit", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{method}{
    string representing the baseline correction method. It can be one of these methods:
    \itemize{
      \item "als" Asymmetric Least Squares, baseline correction by 2nd derivative constrained weighted regression
      \item "fillPeaks" An iterative algorithm using suppression of baseline by means in local windows
      \item "irls" Iterative Restricted Least Squares, an algorithm with primary smoothing and repeated baseline suppressions and regressions with 2nd derivative constraint
      \item "lowpass" Low-pass filter, an algorithm for removing baselines based on Fast Fourier Transform filtering
      \item "medianWindow" an implementation and extension of Mark S. Friedrichs' model-free algorithm
      \item "modpolyfit" Modified polynomial fitting, an implementation of Chad A. Lieber and Anita Mahadevan-Jansen's algorithm for polynomial fitting
      \item "peakDetection" A translation from Kevin R. Coombes et al.'s MATLAB code for detecting peaks and removing baselines
      \item "rfbaseline" Robust Baseline Estimation, Wrapper for Andreas F. Ruckstuhl, Matthew P. Jacobson, Robert W. Field, James A. Dodd's algorithm based on LOWESS and weighted regression
      \item "rollingBall" Ideas from Rolling Ball algorithm for X-ray spectra by M. A. Kneen and H. J. Annegarn. Variable window width has been left out
    }
}
  \item{\dots}{
  Additional parameters of the baseline correction method.
}
}

\value{
  Returns the dataset with the data's baseline corrected.
}

\examples{
  ## Example of baseline correction
  data(cassavaPPD)
  dataset.corrected = baseline_correction(cassavaPPD, method = "modpolyfit")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ baseline }
\keyword{ correction }% __ONLY ONE__ keyword per line
