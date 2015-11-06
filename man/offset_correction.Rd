\name{offset_correction}
\alias{offset_correction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Offset correction
}
\description{
  Perform offset correction on the data.
}
\usage{
offset_correction(dataset)
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
  ## Example of offset correction
  data(cassavaPPD)
  cassava.corrected = offset_correction(cassavaPPD)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ offset }
\keyword{ correction }% __ONLY ONE__ keyword per line
