\name{get_samples_names_dx}
\alias{get_samples_names_dx}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Get sample's names from DX files
}
\description{
  Function to get the names of the DX files from a folder.
}
\usage{
get_samples_names_dx(foldername)
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
  ## Example of getting DX sample's names
  dx.sample.names = get_samples_names_dx("data")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ dx }
\keyword{ sample }% __ONLY ONE__ keyword per line
