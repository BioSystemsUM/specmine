\name{scaling_samples}
\alias{scaling_samples}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Scale data matrix
}
\description{
  Performs scaling according to a method.
}
\usage{
scaling_samples(datamat, method = "auto")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datamat}{
   data matrix.
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
 Returns the data matrix scaled.
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ scaling }
\keyword{ data }% __ONLY ONE__ keyword per line
