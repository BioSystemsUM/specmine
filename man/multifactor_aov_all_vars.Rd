\name{multifactor_aov_all_vars}
\alias{multifactor_aov_all_vars}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Multifactor ANOVA
}
\description{
  Perform multi-factor ANOVA on all variables with the selected metadata variables.
}
\usage{
multifactor_aov_all_vars(dataset, metadata.vars, combination)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{metadata.vars}{
    metadata variables to use in ANOVA.
}
  \item{combination}{
    a formula specifying the model.
}
}
\value{
  List where each element is the multifactor anova result of a variable on the dataset.
}


\examples{
  ## Example of multifactor ANOVA on all variables 
  data(propolis)
  propolis = missingvalues_imputation(propolis, "value", value = 0.00005)
  m.aov.results = multifactor_aov_all_vars(propolis, 
		  c("seasons","agroregions"), "seasons*agroregions")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ anova }
\keyword{ multifactor }% __ONLY ONE__ keyword per line
