\name{remove_data_variables}
\alias{remove_data_variables}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Remove data variables
}
\description{
  Remove data variables from the dataset.
}
\usage{
remove_data_variables(dataset, variables.to.remove, 
by.index = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{variables.to.remove}{
    vector with the variables' names to remove.
}
  \item{by.index}{
    if TRUE, the values of the variables.to.remove argument are indexes.  
}
}

\value{
  Returns the dataset with the specified data variables removed.
}

\examples{
  ## Example of removing data variables
  library(specmine.datasets)
  data(cachexia)
  dataset = remove_data_variables(cachexia, c("Creatine","Serine"), 
	    by.index = FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ remove }
\keyword{ variable }% __ONLY ONE__ keyword per line
