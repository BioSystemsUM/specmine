\name{log_transform}
\alias{log_transform}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Logarithmic transformation.
}
\description{
  Performs logarithmic transformation on the data matrix.
}
\usage{
log_transform(datamat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datamat}{
    data matrix.
}
}

\value{
 Returns the data matrix with the logarithmic transformation applied.
}

\examples{
  ## Example of logarithmic transformation
  library(specmine.datasets)
  data(propolis)
  propolis_proc = missingvalues_imputation(propolis)
  datamat.log = log_transform(propolis_proc$data)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ log }
\keyword{ transformation }% __ONLY ONE__ keyword per line
