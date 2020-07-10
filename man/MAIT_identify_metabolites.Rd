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
   xcmsSet object that can be passed. Stored in dataset$xSet.
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

\details{
  After runing the MAIT_identify_metabolites function, the results table can be aaccessed by:
  
  mait.metab.table = mait.metabolites@FeatureInfo@metaboliteTable
  
  where 'mait.metabolites' is the result obtained from runing MAIT_identify_metabolites.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ MAIT }
\keyword{ metabolite }% __ONLY ONE__ keyword per line
