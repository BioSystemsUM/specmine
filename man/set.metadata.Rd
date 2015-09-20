\name{set.metadata}
\alias{set.metadata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Set new metadata
}
\description{
Updates the dataset's metadata with a new one.
}
\usage{
set.metadata(dataset, new.metadata)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
  list representing the dataset from a metabolomics experiment.
}
  \item{new.metadata}{
  matrix or dataframe with the new metadata.
}
}

\value{
  Returns the dataset with the updated metadata.
}

\examples{
\dontrun{
  ## Example of setting a new metadata to the dataset
  dataset = set.metadata(dataset, new.metadata)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~metadata }
\keyword{ ~dataset }% __ONLY ONE__ keyword per line
