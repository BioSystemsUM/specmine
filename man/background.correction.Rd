\name{background.correction}
\alias{background.correction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Background correction
}
\description{
  Perform background correction on the spectra.
}
\usage{
background.correction(dataset)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
}

\value{
  Returns the dataset with background correction performed on the data.
}

\examples{
\dontrun{
  ## Example of background correction
  dataset.corrected = background.correction(dataset)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~background }
\keyword{ ~correction }% __ONLY ONE__ keyword per line
