\name{set_x_label}
\alias{set_x_label}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Set x-label
}
\description{
  Set a new x-label to the dataset.
}
\usage{
set_x_label(dataset, new.x.label)
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
  ## Example of setting a new x-label to the dataset
  data(cachexia)
  cachexia = set_x_label(cachexia, "new x-label")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ dataset }
\keyword{ label }% __ONLY ONE__ keyword per line
