\name{get_x_values_as_num}
\alias{get_x_values_as_num}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Get x-axis values as numbers
}
\description{
  Get the x-axis values from the dataset as numbers.
}
\usage{
get_x_values_as_num(dataset)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
}

\value{
  Returns a numeric vector with the x-axis values, if the variable labels are not all numeric then an error message is shown.
}

\examples{
  ## Example of getting the x-axis values as numbers
  data(cassavaPPD)
  xvalues.numeric = get_x_values_as_num(cassavaPPD)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ dataset }
\keyword{ xaxis }% __ONLY ONE__ keyword per line
