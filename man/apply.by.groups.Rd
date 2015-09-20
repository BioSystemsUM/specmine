\name{apply.by.groups}
\alias{apply.by.groups}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Apply by groups
}
\description{
 Apply a function to samples from a metadata's variable.
}
\usage{
apply.by.groups(dataset, metadata.var, fn.to.apply = "mean", variables = NULL, variable.bounds = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{metadata.var}{
    name of the metadata's variable.
}
  \item{fn.to.apply}{
    function to apply (e.g. mean, max, min).
}
  \item{variables}{
    allows to define which variables to calculate the stats (if numbers, indexes are assumed).
}
  \item{variable.bounds}{
    allow to define an interval of variables (if numeric).
}
}

\value{
   Returns a vector with the variables and the value of the applied function on the metadata's groups.
}

\examples{
   \dontrun{
     ## Example of applying a function to groups
     apply.groups.result = apply.by.groups(dataset, "metadatavar1", mean)
   }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~groups }
\keyword{ ~apply }% __ONLY ONE__ keyword per line
