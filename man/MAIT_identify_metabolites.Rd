\name{MAIT_identify_metabolites}
\alias{MAIT_identify_metabolites}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  MAIT metabolite identification
}
\description{
 Performs metabolite identification using MAIT.
}
\usage{
MAIT_identify_metabolites(dataset, metadata.variable, 
xSet = NULL, data.folder = NULL, features = NULL, 
mass.tolerance = 0.5)
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
  data(spinalCord)
  library(MAIT)
  mait.metabolites = MAIT_identify_metabolites(spinalCord, "type", 
			features = "all", data.folder = "data",
			xSet = spinalCord$xSet)
  mait.metab.table = mait.metabolites@FeatureInfo@metaboliteTable
  mait.metab.table[which(mait.metab.table$Name != "Unknown"), 
  c(1,3,6)]
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ MAIT }
\keyword{ metabolite }% __ONLY ONE__ keyword per line
