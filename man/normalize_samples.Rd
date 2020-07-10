\name{normalize_samples}
\alias{normalize_samples}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Normalize samples
}
\description{
  Normalize the data from a datamatrix with a specific method.
}
\usage{
normalize_samples(datamat, method, ref = NULL, constant = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datamat}{
    data matrix.
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
  Returns the data matrix normalized.
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ normalization }
\keyword{ sample }% __ONLY ONE__ keyword per line
