\name{convert_to_hyperspec}
\alias{convert_to_hyperspec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Convert to hyperspec
}
\description{
  Convert a dataset to an hyperspec object.
}
\usage{
convert_to_hyperspec(dataset)
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
  ## Example of converting a dataset to an hyperspec object
  data(cassavaPPD)
  hyperspec.cassava = convert_to_hyperspec(cassavaPPD) 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hyperspec }
\keyword{ dataset }% __ONLY ONE__ keyword per line
