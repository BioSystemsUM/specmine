\name{first.derivative}
\alias{first.derivative}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
First derivative
}
\description{
  Calculates the first derivative of the data.
}
\usage{
first.derivative(dataset)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
}

\value{
  Return the dataset with the first derivative of the data calculated.
}

\examples{
\dontrun{
  ## Example of calculate the first derivative
  dataset.derivative = first.derivative(dataset)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~derivative }
\keyword{ ~dataset }% __ONLY ONE__ keyword per line
