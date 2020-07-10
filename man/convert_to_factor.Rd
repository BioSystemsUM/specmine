\name{convert_to_factor}
\alias{convert_to_factor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Convert metadata to factor
}
\description{
  Convert a metadata's variable to factor.
}
\usage{
convert_to_factor(dataset, metadata.var)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
  list representing the dataset from a metabolomics experiment.
}
  \item{metadata.var}{
  name of the metadata's variable.
}
}

\value{
  Returns the dataset with the metadata's variable converted to factor.
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ metadata }
\keyword{ factor }% __ONLY ONE__ keyword per line
