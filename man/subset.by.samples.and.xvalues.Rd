\name{subset.by.samples.and.xvalues}
\alias{subset.by.samples.and.xvalues}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Subset by samples and x-values
}
\description{
 Gets a subset of specific samples and x-values.
}
\usage{
subset.by.samples.and.xvalues(dataset, samples, variables = NULL, by.index = F, variable.bounds = NULL, rebuild.factors = T)
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
\dontrun{
  ## Example of subsetting samples and x-values
  subset = subset.by.samples.and.xvalues(dataset, c("sample1","sample2"), variables = c("var1","var2"))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~subset }
\keyword{ ~sample }% __ONLY ONE__ keyword per line
\keyword{ ~xvalue }