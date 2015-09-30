\name{read_ms_spectra}
\alias{read_ms_spectra}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Read MS spectra
}
\description{
 Read the data from the MS files and creates the dataset.
}
\usage{
read_ms_spectra(folder.name, type = "undefined", 
filename.meta = NULL, description = "", prof.method = "bin", 
fwhm = 30, bw = 30, intvalue = "into", header.col.meta = TRUE, 
header.row.meta = TRUE, sep.meta = ",")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{folder.name}{
   string containing the path of the data folder.
}
  \item{type}{
   type of the data.
}
  \item{filename.meta}{
   name of the metadata file.
}
  \item{description}{
  a short text describing the dataset.
}
  \item{prof.method}{
  profmethod parameter from xcmsSet function from xcms package.
}
  \item{fwhm}{
  fwhm parameter from xcmsSet function from xcms package. A commonly used value is 30 (seconds) for LC-MS and 4 (seconds) for GC-MS spectra.
}
  \item{bw}{
  bw parameter from group function from xcms package.
}
  \item{intvalue}{
  value parameter from groupval function from xcms package. It can be:
    \itemize{
 	\item \bold{"into"} - integrated area of original (raw) peak
	\item \bold{"intf"} - integrated area of filtered peak.
	\item \bold{"maxo"} - maximum intensity of original (raw) peak.
	\item \bold{"maxf"} - maximum intensity of filtered peak.
    }
}
  \item{header.col.meta}{
  boolean value indicating if the metadata CSV file contains a header column with the name of the metadata variables.
}
  \item{header.row.meta}{
  boolean value indicating if the metadata CSV file contains a header row with the name of the samples.
}
  \item{sep.meta}{
  the separator character of the metadata file.
}
}

\value{
 Returns a dataset from the MS files.
}

\examples{
\dontrun{
  ## Example of reading a dataset from MS files
  dataset = read_ms_spectra("data", type = "nmr-spectra", "metadata.csv", 
	    description = "description of the dataset")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ mass }
\keyword{ spectrometry }% __ONLY ONE__ keyword per line
