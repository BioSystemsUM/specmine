\name{MAIT.identify.metabolites}
\alias{MAIT.identify.metabolites}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  MAIT metabolite identification
}
\description{
 Performs metabolite identification using MAIT.
}
\usage{
MAIT.identify.metabolites(dataset, metadata.variable, xSet = NULL, data.folder = NULL, features = NULL, mass.tolerance = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{metadata.variable}{
    metadata's variable.
}
  \item{xSet}{
   xcmsSet object that can be passed.
}
  \item{data.folder}{
   string indicating the data folder.
}
  \item{features}{
   features that can be used to help to identify the metabolites.
}
  \item{mass.tolerance}{
   mass tolerance.
}
}

\value{
 Returns an object resulted from identifyMetabolites function from MAIT package.
}
\references{
 http://www.bioconductor.org/packages/release/bioc/html/MAIT.html
}

\examples{
  \dontrun{
  ## Example of MAIT metabolite identification
  mait.metabolites = MAIT.identify.metabolites(dataset, "metavar1", features = "all", 
    			data.folder = "data")
 }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~MAIT }
\keyword{ ~metabolite }% __ONLY ONE__ keyword per line
