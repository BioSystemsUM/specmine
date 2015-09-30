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
    metadata variables to use in linear regression.
}
  \item{combination}{
    a formula specifying the model.
}
}

\value{
 Returns a list where each element is the linear regression result of a variable on the dataset.
}

\examples{
\dontrun{
  ## Example of linear regression on all variables 
  linreg.results = linreg_all_vars(metabolomics.dataset, c("metadata1",
		   "metadata2"), "metadata1+metadata2")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ linear }
\keyword{ regression }% __ONLY ONE__ keyword per line
