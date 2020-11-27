\name{group_peaks}
\alias{group_peaks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Group peaks
}
\description{
  Group peaks with peak alignment.
}
\usage{
group_peaks(sample.list, type, method = "own", metadata = NULL, 
samp.classes = 1, description = "", label.x = NULL, 
label.values = NULL, step = 0.03)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sample.list}{
   list containing the sample's data.
}
  \item{type}{
   type of the data.
}
  \item{method}{
   method of peak alignment. Can be "own" or "metaboanalyst", which the later is for using the peak alignment used in MetaboAnalyst software.
}
  \item{metadata}{
   data frame containing the metadata.
}
  \item{samp.classes}{
   the metadata's variable to be used in the MetaboAnalyst method.
}
  \item{description}{
   short description of the data.
}
  \item{label.x}{
  the label for the x values.
}
  \item{label.values}{
  the label for the y values.
}
  \item{step}{
   step value for the peak alignment process.
}
}

\value{
 Returns a dataset with the peaks of the data aligned.
}

\examples{
\donttest{
   ## Example of grouping peaks (computationally heavy)
   library(specmine.datasets)
   data(propolisSampleList)
   peaks.ds = group_peaks(propolisSampleList, "nmr-peaks", method = "own", 
	   	metadata = NULL, description = "short description", 
		label.x = "ppm", label.values = "intensity", step = 0.03)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ peak }
\keyword{ alignment }% __ONLY ONE__ keyword per line
