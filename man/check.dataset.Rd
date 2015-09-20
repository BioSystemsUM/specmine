\name{check.dataset}
\alias{check.dataset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Check dataset
}
\description{
  Check if the dataset is valid and if not give the proper error message.
}
\usage{
check.dataset(dataset)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
}

\value{
  Message saying if the dataset is valid or invalid, in the last case also gives the reason.
}


\examples{
\dontrun{
  ## Example checking a dataset
  check.dataset(dataset)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~dataset }
