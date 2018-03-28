\name{detect_nmr_peaks_from_dataset}
\alias{detect_nmr_peaks_from_dataset}
\title{
Detection of the peaks in an NMR spectra dataset.
}
\description{
This function detects the peaks, that have a minimum intensity of \code{baseline_tresh}, and performs alignment of those peaks. 
}
\usage{
detect_nmr_peaks_from_dataset(dataset, baseline_tresh=50000,
                              ap.method="own", ap.samp.classes=1, ap.step=0.03)
}
\arguments{
  \item{dataset}{list representing the dataset from a metabolomics experiment.}
  
  \item{baseline_tresh}{Minimum intensity value that peaks must have. Peaks with intensity smaller than baseline_tresh will not be considered as detected peaks.}
  
  \item{ap.method}{Method to used in the alignment of peaks, after they are identified. Can be "own" or "metaboanalyst", which the later is for using the peak alignment used in MetaboAnalyst software.}
  
  \item{ap.samp.classes}{the metadata's variable to be used in the MetaboAnalyst method.}
  
  \item{ap.step}{step value for the peak alignment process}
}
\value{
Returns a dataset with the peaks detected and aligned.
}