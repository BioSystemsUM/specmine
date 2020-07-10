\name{read_dataset_csv}
\alias{read_dataset_csv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Read dataset from CSV
}
\description{
  Reads the data from a CSV file and creates the dataset.
}
\usage{
read_dataset_csv(filename.data, filename.meta = NULL, 
type = "undefined", description = "", label.x = NULL, 
label.values = NULL, sample.names = NULL, format = "row", 
header.col = TRUE, header.row = TRUE, sep = ",", 
header.col.meta = TRUE, header.row.meta = TRUE, sep.meta = ",")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename.data}{
  name of the data file.
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
  \item{label.x}{
  the label for the x values.
}
  \item{label.values}{
  the label for the y values.
}
  \item{sample.names}{
  the names of the samples.
}
  \item{format}{
  format which the data are in the CSV file. It can be "row" if the samples are in the rows or "col" if the samples are in the columns.
}
  \item{header.col}{
  boolean value indicating if the CSV contains a header column with the names of the  samples or variables.
}
  \item{header.row}{
  boolean value indicating if the CSV contains a header row with the names of the samples or variables.
}
  \item{sep}{
  the separator character.
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
  Returns the dataset from the CSV file.
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ dataset }
\keyword{ csv }% __ONLY ONE__ keyword per line
