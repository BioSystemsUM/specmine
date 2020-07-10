\name{subset_by_samples_and_xvalues}
\alias{subset_by_samples_and_xvalues}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Subset by samples and x-values
}
\description{
 Gets a subset of specific samples and x-values.
}
\usage{
subset_by_samples_and_xvalues(dataset, samples, variables = NULL, 
by.index = FALSE, variable.bounds = NULL, rebuild.factors = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{samples}{
    vector with indexes or names of the samples to select
}
  \item{variables}{
    vector with the desired variables to get from the dataset.
}
  \item{by.index}{
    if TRUE, the values of the variables argument are indexes.  
}
  \item{variable.bounds}{
    variable bounds used if by.index is FALSE and variables are NULL.
}
  \item{rebuild.factors}{
    if TRUE the metadata factors are rebuilded.
}
}

\value{
  Returns the dataset with the selected samples and x-values.
}


\examples{
  ## Example of subsetting samples and x-values
  data(cachexia)
  subset = subset_by_samples_and_xvalues(cachexia, c("PIF_178","NETL_022_V1"), 
	   variables = c("Creatine","Serine"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ subset }
\keyword{ sample }% __ONLY ONE__ keyword per line
\keyword{ xvalue }
