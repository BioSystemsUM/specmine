\name{peak_detection2d}
\alias{peak_detection2d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Detection of the peaks in an 2D NMR spectra dataset.
}
\description{
  This function detects the peaks across samples, reducing the dimensionality of the 2D spectra. If \code{baseline_thresh} is not provided it will be calculated a signal-to-noise ratio (SNR) for each spectra that will serve as threshold. 
}
\usage{
peak_detection2d(specmine_2d_dataset, baseline_thresh=NULL, noiseFilt=0, negatives=F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{specmine_2d_dataset}{
  List representing a 2D dataset from a 2D metabolomics experiment.
  }
  
  \item{baseline_thresh}{
  Minimum intensity value that peaks must have. Peaks with intensity smaller than \code{baseline_thresh} will not be considered as detected peaks.
  }
  
  \item{noiseFilt}{
  Integer argument that applies a noise filter when searching for peaks. Can be one of the following:
    \itemize{
      \item{0}{
      Does not apply a noise filter.
      }
      \item{1}{
      Applies a mild filter (adjacent points in the direct dimension must be above the noise threshold).
      }
      \item{2}{
      Applies a strong filter (all adjacent points must be above the noise threshold).
      }
    }
  Defaults to 0.
  }
  
  \item{negatives}{
  Boolean value that decides if negative ppm values should be considered or not. Default to FALSE.
  }
  
}

\value{
  Returns a 1D specmine dataset where the rows are combinations of ppms (indirect dimension x direct dimension) and the columns are the samples. The combinations of ppms represent peaks detected.
}


