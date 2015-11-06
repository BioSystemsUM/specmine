\name{background_correction}
\alias{background_correction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Background correction
}
\description{
  Perform background correction on the spectra.
}
\usage{
background_correction(dataset)
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
  ## Example of background correction
  library(hyperSpec)
  data(flu)
  flu.converted = convert_from_hyperspec(flu)
  flu.corrected = background_correction(flu.converted)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ background }
\keyword{ correction }% __ONLY ONE__ keyword per line
