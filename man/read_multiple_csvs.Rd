\name{read_multiple_csvs}
\alias{read_multiple_csvs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read multiple CSVs
}
\description{
  Reads multiple CSVs, each one with a sample.
}
\usage{
read_multiple_csvs(filenames, ext = ".csv", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filenames}{
    list of file names of the files to read.
}
  \item{ext}{
    extension name.
}
  \item{\dots}{
    additional parameters to read.csv function.
}
}

\value{
  returns a list of dataframes.
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ read }
\keyword{ csv }% __ONLY ONE__ keyword per line
