\name{get_metadata_value}
\alias{get_metadata_value}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Get metadata value
}
\description{
  Get the metadata value
}
\usage{
get_metadata_value(dataset, variable, sample)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
list representing the dataset from a metabolomics experiment.
}
  \item{variable}{
index or name of the metadata variable.
}
  \item{sample}{
index or name of the sample.
}
}

\value{
Return the corresponding metadata value of the sample.
}

\examples{
\dontrun{
  ## Example of getting a metadata value
  metadata.value = get_metadata_value(dataset, "var1", "sample2")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ metadata }
\keyword{ value }% __ONLY ONE__ keyword per line
