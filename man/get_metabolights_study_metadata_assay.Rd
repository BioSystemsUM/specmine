\name{get_metabolights_study_metadata_assay}
\alias{get_metabolights_study_metadata_assay}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Download metadata file from an assay of MetaboLights study
}
\description{
Downloads the metadata file from the assay specified in 'assay' of the MetaboLights study ('studyID').
}
\usage{
get_metabolights_study_metadata_assay(studyID, assay, directory)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{studyID}{ID of the metabolights study to download.}
  \item{assay}{Number of the assay.}
  \item{directory}{Directory where to download the data.}
}
\details{
This function should be used together with \code{\link{get_metabolights_study_files_assay}}. See example below.

The metadata file is a csv file with the metadata information on each sample of the study's assay.
}
\references{
MetaboLights database: https://www.ebi.ac.uk/metabolights/
}
