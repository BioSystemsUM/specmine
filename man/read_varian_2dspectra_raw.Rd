\name{read_varian_2dspectra_raw}
\alias{read_varian_2dspectra_raw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Function that reads raw 2D spectra (intensity over time spectra) from the varian format and processes them to ppm spectra.
}
\description{
  This function read raw 2D spectra (i.e. intensity over time spectra) from the varian format and processess them to intensity over ppm spectra. For this function to work, in each spectrum directory should be present a fid and procpar files.
  
  Python 3 with module nmrglue must be installed.
}
\usage{
read_varian_2dspectra_raw(varian_spectra_directory,
                          metadata_file=NULL, m.header_col=T, m.header_row=T, m.sep=",",
                          samples.names=NULL, zero_filling=T, apodization=T, zipped=T,
                          description="", label.x="ppm", 
                          label.y="ppm", label.values="intensity")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{varian_spectra_directory}{
  Path of the directory with all the directories of the varian 2D spectra.
  }
  
  \item{metadata_file}{
  Path of the metadata file.
  }
  
  \item{m.header_col}{
  Boolean value indicating if the metadata CSV file contains a header column with the name of the metadata variables.
  }
  
  \item{m.header_row}{
  Boolean value indicating if the metadata CSV file contains a header row with the name of the samples.
  }
  
  \item{m.sep}{
  The separator character of the metadata file.
  }
  
  \item{samples.names}{
  CSV file where the first column represents the samples names and in the second column the names of the spectra directories to which they correspond. If NULL, it will be considered that the directories names are the samples names (it has to be the same names that appear in the metadata file).
  }
  
  \item{zero_filling}{
  Boolean value indicating whether zero-filling should be performed or not when processing the fid of the 2D spectra. Defaults to TRUE.
  }
  
  \item{apodization}{
  Boolean value indicating whether \bold{exponential} apodization should be performed or not when processing the fid of the 2D spectra. Defaults to TRUE. 
  }
  
  \item{zipped}{
  Boolean value indicating if the spectra directories are zipped or not. Defaults to TRUE.
  }
  
  \item{description}{
  A short text describing the dataset.
  }
  
  \item{label.x}{
  The label for the x values (indirect dimension).
  }
  
  \item{label.y}{
  The label for the y values (direct dimension).
  }
  
  \item{label.values}{
  The label for the pair'(x,y) values.
  }
}

\value{
  Returns a list representing a 2D dataset for specmine:
  \item{data}{
  A list of matrices where each matrix matches one 2D spectra.
  }
  \item{type}{
  The type of the data in the dataset.
  }
  \item{description}{
  A short text describing the dataset.
  }
  \item{metadata}{
  A dataframe with the metadata variables.
  }
  \item{F1_ppm}{
  The ppm values regarding indirect dimension.
  }
  \item{F2_ppm}{
  The ppm values regarding direct dimension.
  }
  \item{labels}{
  A list of vectors for the x, y and pairs'(x,y) values.
  }
}

\section{Warning}{
You must not call this function unless you have Python (>=3.5.2) installed in your machine and the module nmrglue.
}
