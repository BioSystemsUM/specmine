\name{read_data_spc}
\alias{read_data_spc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Read data from SPC files
}
\description{
  Reads the data from the SPC files.
}
\usage{
read_data_spc(foldername, nosubhdr = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{foldername}{
   string with the path of the data folder.
}
  \item{nosubhdr}{
   if TRUE, it won't read the subheader.
}
}

\value{
  Returns a list with the samples and the respective data points.
}

\examples{
\dontrun{
  ## Example of reading a dataset from SPC files
  sample.list = read_data_spc("data")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ spc }
\keyword{ dataset }% __ONLY ONE__ keyword per line
