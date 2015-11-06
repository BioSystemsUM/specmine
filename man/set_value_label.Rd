\name{set_value_label}
\alias{set_value_label}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Set value label
}
\description{
  Set a new value label for the dataset.
}
\usage{
set_value_label(dataset, new.val.label)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
  list representing the dataset from a metabolomics experiment.
}
  \item{new.val.label}{
  string with the new value label.
}
}

\value{
  Returns the dataset with the updated value label.
}

\examples{
  ## Example of setting a new value label to the dataset
  data(cachexia)
  cachexia = set_value_label(cachexia, "new value label")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ label }
\keyword{ dataset }% __ONLY ONE__ keyword per line
