\name{scaling}
\alias{scaling}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Scale dataset
}
\description{
  Performs scaling according to a method.
}
\usage{
scaling(dataset, method = "auto")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
      list representing the dataset from a metabolomics experiment.
}
  \item{method}{
    string specifying the scaling method. The possible values are:
      \itemize{
        \item \bold{"auto"} auto scaling.
        \item \bold{"range"} range scaling.
        \item \bold{"pareto"} pareto scaling.
        \item \bold{"tointerval"} scaling to an interval.
      }
}
}

\value{
 Returns the dataset scaled.
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ scaling }
\keyword{ dataset }% __ONLY ONE__ keyword per line
