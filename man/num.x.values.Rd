\name{num.x.values}
\alias{num.x.values}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Get number of x values
}
\description{
  Get the number of x-axis values.
}
\usage{
num.x.values(dataset)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
}

\value{
  Returns an integer with the number of x-axis values.
}

\examples{
\dontrun{
  ## Example of getting the number of x-axis values
  number.x.values = num.x.values(dataset)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~dataset }
\keyword{ ~xaxis }% __ONLY ONE__ keyword per line
