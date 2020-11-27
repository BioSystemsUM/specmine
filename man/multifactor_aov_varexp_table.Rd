\name{multifactor_aov_varexp_table}
\alias{multifactor_aov_varexp_table}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Multifactor ANOVA variability explained table
}
\description{
  Gets the variability explained table from the multifactor ANOVA results.
}
\usage{
multifactor_aov_varexp_table(multifactor.aov.results, 
write.file = FALSE, file.out = "multi-anova-varexp.csv")
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
 Returns a data.frame with the variability explained.
}

\examples{
  ## Example of multifactor ANOVA variability explained table
  library(specmine.datasets)
  data(propolis)
  propolis = missingvalues_imputation(propolis, "value", value = 0.00005)
  m.aov.results = multifactor_aov_all_vars(propolis, 
		  c("seasons","agroregions"), "seasons*agroregions")
  m.aov.varepx = multifactor_aov_varexp_table(m.aov.results)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ varexp }
\keyword{ anova }% __ONLY ONE__ keyword per line
