\name{stats_by_variable}
\alias{stats_by_variable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Statistics of variables
}
\description{
  Get a summary of statistics of the variables.
}
\usage{
stats_by_variable(dataset, variables = NULL, 
variable.bounds = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{variables}{
    allows to define which variables to calculate the stats (if numbers, indexes are assumed).
}
  \item{variable.bounds}{
    allow to define an interval of variables (if numeric).
}
}

\value{
  Returns a vector with the a summary of statistics of the variables. 
}

\examples{
     ## Example of getting stats of variables
	 library(specmine.datasets)
     data(cachexia)
     variable.stats.result = stats_by_variable(cachexia)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ statistic }
\keyword{ variable }% __ONLY ONE__ keyword per line
