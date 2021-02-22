\name{sum_2d_dataset}
\alias{sum_2d_dataset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  2D Dataset summary
}
\description{
  Returns a summary with its main features
}
\usage{
sum_2d_dataset(dataset_2d, stats = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset_2d}{
  List representing the 2D dataset from a 2D metabolomics experiment.
  }
  
  \item{stats}{
  If TRUE prints some global statistics for each 2D spectra.
  }
}

\value{
  Returns the summary of the 2D dataset containing:
  \itemize{
    \item Description
    \item Type of data
    \item Number of samples
    \item Number of data points
    \item Number of metadata variables if metadata not null
    \item Labels of x axis, y axis and pair'(x,y) values if labels not null
  }
  
  If the parameter 'stats' is TRUE then returns also:
  \itemize{
    \item Number of missing values in each spectra
    \item Mean of data values in each spectra
    \item Median of data values in each spectra
    \item Standard deviation in each spectra
    \item Range of values in each spectra
    \item Quantiles for each spectra
  }
}



% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ 2d-dataset }
\keyword{ summary }
