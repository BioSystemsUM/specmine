\name{replace_metadata_value}
\alias{replace_metadata_value}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Replace metadata's value
}
\description{
  Replace a metadata's variable value of a sample.
}
\usage{
replace_metadata_value(dataset, variable, sample, new.value)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
  list representing the dataset from a metabolomics experiment.
}
  \item{variable}{
  metadata's variable.
}
  \item{sample}{
  name of the sample.
}
  \item{new.value}{
  new value of the metadata's variable.
}
}

\value{
  Returns the dataset with the metadata updated.
}

\examples{
  ## Example of replacing metadata's variable value of a sample
  data(cachexia)
  dataset = replace_metadata_value(cachexia, "Muscle.loss", "PIF_178", 
	    "control")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ metadata }
\keyword{ dataset }% __ONLY ONE__ keyword per line
