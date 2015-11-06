\name{replace_data_value}
\alias{replace_data_value}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Replace data value
}
\description{
  Replace a data value for a new value on the dataset.
}
\usage{
replace_data_value(dataset, x.axis.val, sample, new.value, 
by.index = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
  list representing the dataset from a metabolomics experiment.
}
  \item{x.axis.val}{
  variable index or name.
}
  \item{sample}{
  sample name.
}
  \item{new.value}{
  new value to replace the old value.
}
  \item{by.index}{
  boolean value to indicate if the x.axis.val is an index or not.
}
}

\value{
  Returns the dataset with the data value replaced.
}

\examples{
  ## Example of replacing a data value from the dataset
  data(cachexia)
  dataset = replace_data_value(cachexia, "Creatine", "PIF_178", 10.3, 
	    by.index = FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ dataset }
\keyword{ value }% __ONLY ONE__ keyword per line
