\name{get_value_label}
\alias{get_value_label}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Get value label
}
\description{
  Get the value label from the dataset
}
\usage{
get_value_label(dataset)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
}

\value{
  Returns a string with the value label.
}


\examples{
  ## Example of getting the value label
  data(cassavaPPD)
  value.label = get_value_label(cassavaPPD)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ label }
\keyword{ dataset }% __ONLY ONE__ keyword per line
