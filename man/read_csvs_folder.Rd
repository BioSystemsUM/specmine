\name{read_csvs_folder}
\alias{read_csvs_folder}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Read CSVs from folder
}
\description{
  Reads multiple CSV files in a given folder.
}
\usage{
read_csvs_folder(foldername, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{foldername}{
    string with the name of the folder.
}
  \item{\dots}{
    additional parameters to read.csv function.
}
}

\value{
  Returns a list of data frames.
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ read }
\keyword{ csv }% __ONLY ONE__ keyword per line
