\name{remove.data}
\alias{remove.data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Remove data
}
\description{
 Remove data from the dataset.
}
\usage{
remove.data(dataset, data.to.remove, type = "sample", by.index = F, rebuild.factors = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{data.to.remove}{
    vector with the elements' names to remove
}
  \item{type}{
    type of the element to remove. It can be:
    \itemize{
    	\item \bold{"sample"} to remove samples.
        \item \bold{"data"} to remove variables.
        \item \bold{"metadata"} to remove metadata's variables.
    }
}
  \item{by.index}{
    if TRUE, the values of the data.to.remove argument are indexes in case the type is "data".  
}
  \item{rebuild.factors}{
    if TRUE, rebuilds the factors from metadata.
}
}

\value{
  Returns the dataset with the specified data removed.
}


\examples{
\dontrun{
  ## Example of removing data
  dataset = remove.data(dataset, c("var1","var2"), type = "data", by.index = F)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~remove }
\keyword{ ~data }% __ONLY ONE__ keyword per line
