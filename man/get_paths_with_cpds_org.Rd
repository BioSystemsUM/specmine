\name{get_paths_with_cpds_org}
\alias{get_paths_with_cpds_org}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get only the paths of the organism that contain one or more of the given compounds.
}
\description{
Gives only the metabolic paths of the mentioned organism that contain one or more of the given compounds.
}
\usage{
get_paths_with_cpds_org(organism_code, compounds, full.result=T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{organism_code}{Organism code. The correct code for an organism can be consulted using function \code{\link{get_OrganismsCodes}}.}
  \item{compounds}{Named vector with kegg codes of compounds and respective names. This vector can be obtained by using the function \code{\link{get_cpd_names}} or the function \code{\link{convert_hmdb_to_kegg}}.}
  \item{full.result}{If the full result is to be given. Defaults to TRUE.}
}
\value{
Data frame.

If full result is chosen, the data frame contains information on the pathways of the organism that contains one or more of the given compounds and, for each pathway, the kegg codes (and their names) of the compounds given that are present in that path. The ratio between the number of compounds given compounds present in each pathway and the total number of compounds in each pathway is also given, full result or not.

If full result is not wanted, only the pathways will be given.
}
\examples{
#Get human metabolic paths that have one or more of the three following compounds
#keggs=get_cpd_names(c("cpd:C00033", "cpd:C00147", "gl:G13099"))
#paths_org_cpds=get_paths_with_cpds_org("hsa", keggs)
#paths_org_cpds
}
