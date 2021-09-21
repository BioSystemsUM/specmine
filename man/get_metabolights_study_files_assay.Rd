\name{get_metabolights_study_files_assay}
\alias{get_metabolights_study_files_assay}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Download data files from an assay of MetaboLights study
}
\description{
Downloads the data files from the assay specified in 'assay' of the MetaboLights study ('studyID')
}
\usage{
get_metabolights_study_files_assay(studyID, assay, directory)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{studyID}{ID of the metabolights study to download.}
  \item{assay}{Number of the assay.}
  \item{directory}{Directory where to download the data.}
}
\details{
This function should be used together with \code{\link{get_metabolights_study_metadata_assay}}. See example below.

Be aware that the study's files may not be structured in the right way to be readily imported with a specmine read function.

Specmine takes into consideration that the names of the data files/folders correspond to the names of the samples. In some studies, data file names do not correspond to the samples' names in the metadata. To overcome this, we create a file called 'samples_files.csv' matching the sample name to the respective data file/zipped folder.

In some cases, one downloaded zipped data folder may contain more than one sample / replicates, but metabolights information only associates the overall folder as one sample. So manual naming of the folder samples and further changing the metadata file (metadata.csv) may be necessary.

Also, some data formats of some metabolights studies are not yet readable by specmine.
}
\references{
MetaboLights database: https://www.ebi.ac.uk/metabolights/
}
