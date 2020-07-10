\name{linreg_rsquared}
\alias{linreg_rsquared}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Linear regression r-squared
}
\description{
 Gets the linear regression r-squared values.
}
\usage{
linreg_rsquared(linreg.results, write.file = FALSE, 
file.out = "linreg-rsquared.csv")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{linreg.results}{
   Linear regression results from linreg.all.vars function.
}
  \item{write.file}{
   boolean value to indicate if a file should be written with the results.
}
  \item{file.out}{
   name of the file.
}
}

\value{
 Returns a data.frame with the r-squared values.
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ rsquared }
\keyword{ regression }% __ONLY ONE__ keyword per line
