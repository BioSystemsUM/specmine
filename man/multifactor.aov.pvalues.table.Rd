\name{multifactor.aov.pvalues.table}
\alias{multifactor.aov.pvalues.table}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Multifactor ANOVA p-values table
}
\description{
 Gets the p-values table from the multifactor ANOVA results.
}
\usage{
multifactor.aov.pvalues.table(multifactor.aov.results, write.file = F, file.out = "multi-anova-pvalues.csv")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{multifactor.aov.results}{
   multifactor anova results.
}
  \item{write.file}{
   boolean value to indicate if a file is written.
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
  ## Example of multifactor ANOVA p-values table
  m.aov.pvalues = multifactor.aov.pvalues.table(multifactor.aov.results)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~multifactor }
\keyword{ ~anova }% __ONLY ONE__ keyword per line
