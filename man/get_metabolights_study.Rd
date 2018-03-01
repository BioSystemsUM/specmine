\name{get_metabolights_study}
\alias{get_metabolights_study}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Download Metabolights study files.
}
\description{
Download data files and metadata files from the specified study form metabolights databaes.
}
\usage{
get_metabolights_study(mtblsID, directory)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mtblsID}{ID of the metabolights study to download.}
  \item{directory}{Directory to where download the data. The data will be stored in "/directory/mtblsID"}
}
\note{
It also downloads a file (\code{samples_files}) with the names of the samples to which each data zip folder corresponds to. This may not be enough when reading the files to a specmine object, as specmine takes into consideration that the names of the data files/folders correspond to the names of the samples and in some cases one downloaded data folder may contain more than one "sample" (replicates, for example) but metabolights information only associates the overall folder as one sample. So manual naming of the folder samples and further changing the metadata file (metadata.csv) may be necessary.

Also, some data formats of some metabolights studies are not yet readable by specmine.
}

\section{Warning}{
You must not call this function unless you have Python (>=3.5.2) installed in your machine and the modules isatools, os, ftplib, glob, logging, pandas, tempfile, shutil, re.
}

\examples{
##get_metabolights_study("MTBLS131", "~/Documents")
}
