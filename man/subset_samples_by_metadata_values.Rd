\name{subset_samples_by_metadata_values}
\alias{subset_samples_by_metadata_values}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Subset samples by metadata values
}
\description{
  Gets a subset of specific samples according to metadata's values from the dataset.
}
\usage{
subset_samples_by_metadata_values(dataset, metadata.varname, 
values)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{metadata.varname}{
    name of the metadata's variable.
}
  \item{values}{
    values of the metadata's variable.
}
}

\value{
  Returns the dataset with the set of samples according to the metadata's values.
}

\examples{
  ## Example of subsetting samples by metadata's values
  data(cassavaPPD)
  subset = subset_samples_by_metadata_values(cassavaPPD, "varieties", 
	   c("BRA","IAC"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ subset }
\keyword{ metadata }% __ONLY ONE__ keyword per line
