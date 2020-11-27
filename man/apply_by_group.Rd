\name{apply_by_group}
\alias{apply_by_group}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Apply by group
}
\description{
  Apply a function to samples from a given metadata's group.
}
\usage{
apply_by_group(dataset, fn.to.apply, metadata.var, var.value)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{fn.to.apply}{
    function to apply (e.g. mean, max, min).
}
  \item{metadata.var}{
    name of the metadata's variable.
}
  \item{var.value}{
    value of the metadata's variable.
}
}

\value{
   Returns a vector with the variables and the value of the applied function.
}

\examples{
     ## Example of applying a function to a group
	 library(specmine.datasets)
     data(cachexia)
     apply.group.result = apply_by_group(cachexia, mean, "Muscle.loss", 
			  "control")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ apply }
\keyword{ group }% __ONLY ONE__ keyword per line
