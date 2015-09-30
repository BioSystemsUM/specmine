\name{read_data_csv}
\alias{read_data_csv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Read CSV data
}
\description{
  Reads the data from the CSV file.
}
\usage{
read_data_csv(filename, format = "row", header.col = TRUE, 
header.row = TRUE, sep = ",")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{
  name of the file with the data.
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
}

\value{
  Returns a numeric matrix with the data.
}
\examples{
\dontrun{
  ## Example of reading a dataset from a CSV file
  data.matrix = read_data_csv("data.csv", format = "row")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ matrix }
\keyword{ data }% __ONLY ONE__ keyword per line
