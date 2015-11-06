\name{set_sample_names}
\alias{set_sample_names}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Set samples names
}
\description{
  Set new samples names to the dataset.
}
\usage{
set_sample_names(dataset, new.sample.names)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
  list representing the dataset from a metabolomics experiment.
}
  \item{new.sample.names}{
  vector with the new samples names.
}
}

\value{
  Returns the dataset with the updated samples names.
}

\examples{
  ## Example of setting a new value label to the dataset
  data(cachexia)
  new.samples.names = as.character(1:77)
  cachexia = set_sample_names(cachexia, new.samples.names)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ sample }
\keyword{ dataset }% __ONLY ONE__ keyword per line
