\name{plot_2d_spectra}
\alias{plot_2d_spectra}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Plot of 2D spectra
}
\description{
  Plot spectra from 2D specmine dataset. 
}
\usage{
plot_2d_spectra(specmine_2d_dataset, title_spectra = "", meta = NULL, spec_samples = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{specmine_2d_dataset}{
  List representing the 2D dataset from a 2D metabolomics experiment.
  }
  
  \item{title_spectra}{
  The title of the plot.
  }
  
  \item{meta}{
  String indicating the metadata's variable. 
  }
  
  \item{spec_samples}{
  Vector with samples' names or numbers, if NULL only four samples are plotted. The two with highest and lowest signal-to-noise ratio (SNR).
  }
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ 2d-spectra }
\keyword{ plot }


