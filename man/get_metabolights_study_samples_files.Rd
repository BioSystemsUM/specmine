\name{get_metabolights_study_samples_files}
\alias{get_metabolights_study_samples_files}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get list of files from an assay of the MetaboLights study and saves it in a csv file.
}
\description{
Get list of files from an assay of the MetaboLights study and saves it in a csv file.
}
\usage{
get_metabolights_study_samples_files(studyID, assay, directory)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{studyID}{ID of the metabolights study}
  \item{assay}{Number of the assay}
  \item{directory}{Directory path where the file will be saved.}
}
\references{
MetaboLights database: https://www.ebi.ac.uk/metabolights/
}

\examples{
\donttest{
 get_metabolights_study_samples_files('MTBLS346',1, tempdir())
 }
}
