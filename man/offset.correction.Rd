\name{offset.correction}
\alias{offset.correction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Offset correction
}
\description{
  Perform offset correction on the data.
}
\usage{
offset.correction(dataset)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
}

\value{
  Returns the dataset 
}

\examples{
\dontrun{
  ## Example of offset correction
  dataset.corrected = offset.correction(dataset)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~offset }
\keyword{ ~correction }% __ONLY ONE__ keyword per line
