\name{normalize}
\alias{normalize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Normalize data
}
\description{
  Normalize the data from the dataset with a specific method.
}
\usage{
normalize(dataset, method, ref = NULL, constant = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
      list representing the dataset from a metabolomics experiment.
}
  \item{method}{
    string specifying the normalization method. The possible values are:
      \itemize{
        \item \bold{"sum"} normalization by sum.
        \item \bold{"median"} normalization by median.
        \item \bold{"ref.sample"} normalization by reference sample.
        \item \bold{"ref.feature"} normalization by reference feature.
      }
}
  \item{ref}{
    the reference if method is "ref.sample" or "ref.feature".
}
  \item{constant}{
    the constant value if method is "sum".
}
}

\value{
  Returns the dataset with the data normalized.
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ normalization }
\keyword{ data }% __ONLY ONE__ keyword per line
