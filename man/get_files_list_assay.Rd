\name{get_files_list_per_assay}
\alias{get_files_list_per_assay}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get list of files per assay for MetaboLights study.
}
\description{
Returns a list of the data files in each assay of a MetaboLights study.
}
\usage{
get_files_list_per_assay(studyID)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{studyID}{ID of the metabolights study}
}
\value{
A list with one or more item. Each item corresponds to an assay of the MetaboLights study. Each item contains a data frame with the names of the samples (column 'Samples') and respective file names (column 'Files').
}
\references{
MetaboLights database: https://www.ebi.ac.uk/metabolights/
}
\examples{
\donttest{
get_files_list_per_assay('MTBLS346')
}
}
