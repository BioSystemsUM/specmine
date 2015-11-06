\name{apply_by_variable}
\alias{apply_by_variable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Apply function to variables
}

\description{
  Applies a function to the values of each variable
}
\usage{
apply_by_variable(dataset, fn.to.apply, variables = NULL, 
variable.bounds = NULL, samples = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
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
  \item{samples}{
    if defined restricts the application to a given set of samples.
}
  \item{\dots}{
    additional parameters to apply function.
}
}

\value{
   Returns a vector with the variables and the value of the applied function.
}


\examples{
     ## Example of applying a function to variables
     data(cachexia)
     apply.variables.result = apply_by_variable(cachexia, mean)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ apply }
\keyword{ variable }% __ONLY ONE__ keyword per line
