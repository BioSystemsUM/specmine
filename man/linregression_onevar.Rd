\name{linregression_onevar}
\alias{linregression_onevar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Linear regression on one variable
}
\description{
 Performs linear regression on one variable of the dataset.
}
\usage{
linregression_onevar(dataset, x.val, metadata.vars, combination)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{x.val}{
    the x-value to be tested.
}
  \item{metadata.vars}{
    metadata variables to use in linear regression. For example, c('variable1','variable2').
}
  \item{combination}{
    a formula specifying the model. For example, 'variable1+variable2'.
}
}

\value{
 Returns a summary of the result from the lm function from stats package.	
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ linear }
\keyword{ regression }% __ONLY ONE__ keyword per line
