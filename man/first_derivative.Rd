\name{first_derivative}
\alias{first_derivative}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
First derivative
}
\description{
  Calculates the first derivative of the data.
}
\usage{
first_derivative(dataset)
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
  ## Example of calculate the first derivative
  data(cassavaPPD)
  cassava.derivative = first_derivative(cassavaPPD)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ derivative }
\keyword{ dataset }% __ONLY ONE__ keyword per line
