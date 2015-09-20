\name{convert.to.hyperspec}
\alias{convert.to.hyperspec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Convert to hyperspec
}
\description{
  Convert a dataset to an hyperspec object.
}
\usage{
convert.to.hyperspec(dataset)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
}

\value{
   Returns an hyperspec object representing the dataset converted.
}


\examples{
\dontrun{
  ## Example of converting a dataset to an hyperspec object
  hyperspec.dataset = convert.to.hyperspec(dataset) 
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~hyperspec }
\keyword{ ~dataset }% __ONLY ONE__ keyword per line
