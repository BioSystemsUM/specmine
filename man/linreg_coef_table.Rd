\name{linreg_coef_table}
\alias{linreg_coef_table}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Linear regression coefficient table
}
\description{
 Gets a data.frame with the coefficient values.
}
\usage{
linreg_coef_table(linreg.results, write.file = F, 
file.out = "linreg-coefs.csv")
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
 Returns a data.frame with the coefficient values.
}

\examples{
\dontrun{
  ## Example of linear regression coefficient table
  linreg.coef.tab = linreg_coef_table(linreg.results)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ coefficient }
\keyword{ regression }% __ONLY ONE__ keyword per line
