\name{get_x_label}
\alias{get_x_label}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Get x-axis label
}
\description{
  Get the x-axis label from the dataset.
}
\usage{
get_x_label(dataset)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
}

\value{
  Returns a string with the x-axis label.
}

\examples{
  ## Example of getting the x-axis label
  data(cassavaPPD)
  x.label = get_x_label(propolis)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ label }
\keyword{ xaxis }% __ONLY ONE__ keyword per line
