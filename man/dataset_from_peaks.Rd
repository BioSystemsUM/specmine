\name{dataset_from_peaks}
\alias{dataset_from_peaks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Dataset from peaks
}
\description{
  Converts a peak list to a dataset.
}
\usage{
dataset_from_peaks(sample.list, metadata = NULL, 
description = "", type = "nmr-peaks")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sample.list}{
    list with the peaks from each sample.
}
  \item{metadata}{
    data frame with the associated metadata.
}
  \item{description}{
    string with the description of the dataset.
}
  \item{type}{
    string that represents the type of the data.
  }
}

\value{
  Returns the dataset from the peak list.
}

\examples{
\dontrun{
    ## Example of converting a peak list to a dataset (computationally heavy)
    data(propolisSampleList)
    dataset = dataset_from_peaks(propolisSampleList, metadata = NULL, 
	      description = "some text", type = "nmr-peaks")  
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ peaklist }
\keyword{ dataset }% __ONLY ONE__ keyword per line
