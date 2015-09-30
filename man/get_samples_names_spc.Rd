\name{get_samples_names_spc}
\alias{get_samples_names_spc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Get sample's names from SPC files
}
\description{
  Function to get the names of the SPC files from a folder.
}
\usage{
get_samples_names_spc(foldername)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{foldername}{
     string with the path of the data folder.
}
}

\value{
  Returns a vector with the sample's names.
}

\examples{
\dontrun{
  ## Example of getting SPC sample's names
  spc.sample.names = get_samples_names_spc("data")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ spc }
\keyword{ names }% __ONLY ONE__ keyword per line
