\name{set.x.label}
\alias{set.x.label}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Set x-label
}
\description{
  Set a new x-label to the dataset.
}
\usage{
set.x.label(dataset, new.x.label)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
  list representing the dataset from a metabolomics experiment.
}
  \item{new.x.label}{
  string with the x-label.
}
}

\value{
  Returns the dataset with the updated x-label.
}

\examples{
\dontrun{
  ## Example of setting a new x-label to the dataset
  dataset = set.x.label(dataset, "new x-label")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~dataset }
\keyword{ ~label }% __ONLY ONE__ keyword per line
