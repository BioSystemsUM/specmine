\name{linreg_all_vars}
\alias{linreg_all_vars}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Linear Regression
}
\description{
 Performs linear regression analysis over the dataset with the selected metadata's variables.
}
\usage{
linreg_all_vars(dataset, metadata.vars, combination)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{metadata.vars}{
    metadata variables to use in linear regression. For example, c('variable1','variable2').
}
  \item{combination}{
    a formula specifying the model. For example, 'variable1+variable2'.
}
}

\value{
 Returns a list where each element is the linear regression result of a variable on the dataset.
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ linear }
\keyword{ regression }% __ONLY ONE__ keyword per line
