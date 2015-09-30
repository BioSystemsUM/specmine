\name{linreg_pvalue_table}
\alias{linreg_pvalue_table}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Linear regression p-values table
}
\description{
 Gets the p-values table from the linear regression analysis.
}
\usage{
linreg_pvalue_table(linreg.results, write.file = F, 
file.out = "linreg-pvalues.csv")
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
 Returns a data.frame with the p-values.
}

\examples{
\dontrun{
  ## Example of linear regression p-values table
  linreg.pvalue.tab = linreg_pvalue_table(linreg.results)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ pvalue }
\keyword{ regression }% __ONLY ONE__ keyword per line
