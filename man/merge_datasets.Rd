\name{merge_datasets}
\alias{merge_datasets}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Merge two datasets
}
\description{
  Merges two datasets with the same variables and metadata's variables.
}
\usage{
merge_datasets(dataset1, dataset2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset1}{
  list representing the first dataset from a metabolomics experiment.
}
  \item{dataset2}{
  list representing the second dataset from a metabolomics experiment.
}
}

\value{
  Returns one dataset with the data from the two datasets merged.
}

\examples{
\dontrun{
  ## Example of merging two datasets
  dataset = merge_datasets(dataset1, dataset2)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ dataset }
\keyword{ merge }% __ONLY ONE__ keyword per line
