\name{read_dataset_spc}
\alias{read_dataset_spc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Read dataset from SPC files
}
\description{
    Reads the data from the SPC files and creates the dataset.
}
\usage{
read_dataset_spc(folder.data, filename.meta = NULL, 
type = "undefined", description = "", nosubhdr = FALSE, 
label.x = NULL, label.values = NULL, header.col.meta = TRUE, 
header.row.meta = TRUE, sep.meta = ",")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{folder.data}{
      string containing the path of the data folder.
}
  \item{filename.meta}{
   name of the metadata file.
}
  \item{type}{
   type of the data.
}
  \item{description}{
  a short text describing the dataset.
}
  \item{nosubhdr}{
  if TRUE, it won't read the subheader.
}
  \item{label.x}{
  the label for the x values.
}
  \item{label.values}{
  the label for the y values.
}
  \item{header.col.meta}{
  boolean value indicating if the metadata CSV file contains a header column with the name of the metadata variables.
}
  \item{header.row.meta}{
  boolean value indicating if the metadata CSV file contains a header row with the name of the samples.
}
  \item{sep.meta}{
  the separator character of the metadata file.
}
}

\value{
  Returns the dataset from the SPC files.
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ dataset }
\keyword{ spc }% __ONLY ONE__ keyword per line
