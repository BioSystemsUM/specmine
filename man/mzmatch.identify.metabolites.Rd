\name{mzmatch.identify.metabolites}
\alias{mzmatch.identify.metabolites}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 mzMatch metabolite identification
}
\description{
 Performs metabolite identification using mzMatch.
}
\usage{
mzmatch.identify.metabolites(ionisation = "detect", data.folder = NULL, xSet = NULL, backend = "Ramp", adducts = "M+H,M+ACN+Na,M+Na,M+K,M+ACN+H", databases = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ionisation}{
   ionisation of the samples. Can be "detect", "positive" or "negative".
}
  \item{data.folder}{
   string indicating the data folder.
}
  \item{xSet}{
   xcmsSet object that can be passed.
}
  \item{backend}{
   backend used to read the files. Can be "Ramp" or "netCDF".
}
  \item{adducts}{
   string with the adducts information.
}
  \item{databases}{
   databases to search for the metabolites.
}
}

\value{
 Returns a list with the final result of the .peakml file read.
}
\references{
  mzmatch.sourceforge.net/
}

\examples{
 \dontrun{
  ## Example of mzmatch metabolite identification
  dbs = dir(paste(find.package("mzmatch.R"), "/dbs", sep=""), full.names = TRUE)
  mzmatch.metabolites = mzmatch.identify.metabolites(ionisation = "positive", "data", 
    			adducts = "M+H,M+ACN+NA", databases = dbs)
 }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~mzmatch }
\keyword{ ~metabolite }% __ONLY ONE__ keyword per line
