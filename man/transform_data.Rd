\name{transform_data}
\alias{transform_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Transform data
}
\description{
  Performs data transformation according to a method.
}
\usage{
transform_data(dataset, method = "log")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
      list representing the dataset from a metabolomics experiment.
}
  \item{method}{
    string specifying the transformation method. The possible values are:
      \itemize{
        \item \bold{"log"} logarithmic transformation.
        \item \bold{"cubicroot"} cubic root transformation.
      }
}
}

\value{
  Returns the dataset with the data transformation applied.
}

\examples{

  ## Example of logarithmic transformation
  library(specmine.datasets)
  data(cachexia)
  dataset.log = transform_data(cachexia, "log")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ log }
\keyword{ cubicroot }% __ONLY ONE__ keyword per line
