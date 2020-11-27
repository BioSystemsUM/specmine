\name{cubic_root_transform}
\alias{cubic_root_transform}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Cubic root transformation
}
\description{
  Performs cubic root transformation on the data matrix.
}
\usage{
cubic_root_transform(datamat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datamat}{
   data matrix.
}
}

\value{
   Returns the data matrix with the cubic root transformation applied.
}

\examples{
  ## Example of cubic root transformation
  library(specmine.datasets)
  data(cachexia)
  datamat.cubic = cubic_root_transform(cachexia$data)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cubicroot }
\keyword{ transformation }% __ONLY ONE__ keyword per line
