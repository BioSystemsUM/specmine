\name{replace.metadata.value}
\alias{replace.metadata.value}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Replace metadata's value
}
\description{
  Replace a metadata's variable value of a sample.
}
\usage{
replace.metadata.value(dataset, variable, sample, new.value)
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
\dontrun{
  ## Example of replacing metadata's variable value of a sample
  dataset = replace.metadata.value(dataset, "metadata.var1", "sample4", "new.value")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~metadata }
\keyword{ ~dataset }% __ONLY ONE__ keyword per line
