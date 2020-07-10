\name{get_metabolights_study}
\alias{get_metabolights_study}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Download MetaboLights study files.
}
\description{
Download data and metadata files for each assay from the specified MetaboLights database study.
}
\usage{
get_metabolights_study(studyID, directory)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{studyID}{ID of the metabolights study to download. For example, 'MTBLS100'}
  \item{directory}{Directory where to download the data.}
}
\note{
Study's files are stored by assay. Data files from assay 1 of the study will be stored in folder '1'.

Be aware that the study's files may not be structured in the right way to be readily imported with a specmine read function.

Specmine takes into consideration that the names of the data files/folders correspond to the names of the samples. In some studies, data file names do not correspond to the samples' names in the metadata. To overcome this, we create a file called 'samples_files.csv' matching the sample name to the respective data file/zipped folder.

In some cases, one downloaded zipped data folder may contain more than one sample / replicates, but metabolights information only associates the overall folder as one sample. So manual naming of the folder samples and further changing the metadata file (metadata.csv) may be necessary.

Also, some data formats of some metabolights studies are not yet readable by specmine.

The metadata file(s) are csv file(s) with the metadata information on each sample. There is one metadata file per assay. Metadata file from assay 1 will be named 'metadata1.csv'.
}
\references{
MetaboLights database: https://www.ebi.ac.uk/metabolights/
}
