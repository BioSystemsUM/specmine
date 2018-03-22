\name{read_Bruker_files}
\alias{read_Bruker_files}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read Bruker processed spectra.
}
\description{
This functions read a directory containing directories where each one corresponds to a bruker spectrum directory.

A CSV file with the names of the samples and to which bruker spectrum directory (directory name) they correspond to should be given, unless directories' names correspond to the samples names.
}
\usage{
read_Bruker_files(bruker_directory, metadata_file=NULL,
                  m.header_col=T, m.header_row=T, m.sep=",",
                  samples.names=NULL, zipped=T,
                  description="", label.x="ppm", label.values="intensity")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bruker_directory}{Path of the directory with all the directories of the bruker spectra.}
  
  \item{metadata_file}{Path of the metadata file.}
  
  \item{m.header_col}{Boolean value indicating if the metadata CSV file contains a header column with the name of the metadata variables.}
  
  \item{m.header_row}{Boolean value indicating if the metadata CSV file contains a header row with the name of the samples.}
  
  \item{m.sep}{The separator character of the metadata file.}
  
  \item{samples.names}{CSV file where the first column represents the samples names and in the second column the names of the spectra directories to which they correspond. If NULL, it will be considered that the directories names are the samples names (it has to be the same names that appear in the metadata file).}
  
  \item{zipped}{Boolean value indicating if the spectra directories are zipped or not.}
  
  \item{description}{A short text describing the dataset.}
  
  \item{label.x}{The label for the x values.}
  
  \item{label.values}{The label for the y values.}
}
\value{
Returns a list representing a dataset for specmine.
}

