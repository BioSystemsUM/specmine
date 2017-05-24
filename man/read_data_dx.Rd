\name{read_data_dx}
\alias{read_data_dx}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Read data from (J)DX files
}
\description{
  Reads the data from the (J)DX files.
}
\usage{
read_data_dx(foldername, debug = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{foldername}{
   string with the path of the data folder.
}
  \item{debug}{
   debug option for readJDX's readJDX function.
}
}

\value{
  Returns a list with the samples and the respective data points.
}

\examples{
\dontrun{
  ## Example of reading a dataset from (J)DX file
s  sample.list = read_data_dx("data")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ jdx }
\keyword{ data }% __ONLY ONE__ keyword per line
