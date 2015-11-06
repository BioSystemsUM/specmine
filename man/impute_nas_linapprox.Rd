\name{impute_nas_linapprox}
\alias{impute_nas_linapprox}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Impute missing values with linear approximation
}
\description{
  Impute missing values with linear approximation.
}
\usage{
impute_nas_linapprox(dataset)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
      list representing the dataset from a metabolomics experiment.
}
}

\value{
Returns the dataset with no missing values.
}

\examples{
  ## Example of NA imputation with linear approximation
  data(propolis)
  dataset = impute_nas_linapprox(propolis)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ missing }
\keyword{ values }% __ONLY ONE__ keyword per line
