\name{set_x_values}
\alias{set_x_values}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Set new x-values
}
\description{
Set new x-values to the dataset
}
\usage{
set_x_values(dataset, new.x.values, new.x.label = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
  list representing the dataset from a metabolomics experiment.
}
  \item{new.x.values}{
  vector with the new x-values.
}
  \item{new.x.label}{
  string with the new x-label (can be NULL).
}
}

\value{
  Returns the dataset with the updated x-values.
}

\examples{
  ## Example of setting new x-values to the dataset
  data(cachexia)
  new.xvalues = 1:63
  cachexia = set_x_values(cachexia, new.xvalues, new.x.label = NULL)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ xvalues }
\keyword{ dataset }% __ONLY ONE__ keyword per line
