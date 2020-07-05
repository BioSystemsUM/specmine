\name{get_x_values_as_text}
\alias{get_x_values_as_text}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Get x-axis values as text
}
\description{
  Get the x-axis values from the dataset as text.
}
\usage{
get_x_values_as_text(dataset)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
}

\value{
  Returns a character vector with the x-axis values.
}


\examples{
  ## Example of getting the x-axis values as text
  xvalues.text = get_x_values_as_text(propolis)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ xaxis }
\keyword{ dataset }% __ONLY ONE__ keyword per line
