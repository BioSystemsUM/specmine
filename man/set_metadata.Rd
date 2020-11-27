\name{set_metadata}
\alias{set_metadata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Set new metadata
}
\description{
Updates the dataset's metadata with a new one.
}
\usage{
set_metadata(dataset, new.metadata)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
  list representing the dataset from a metabolomics experiment.
}
  \item{new.metadata}{
  matrix or dataframe with the new metadata.
}
}

\value{
  Returns the dataset with the updated metadata.
}

\examples{
  ## Example of setting a new metadata to the dataset
  library(specmine.datasets)
  data(cachexia)
  new.metadata = c(rep("meta1", 39), rep("meta2", 38)) 
  new.metadata = data.frame(var_meta = new.metadata)
  rownames(new.metadata) = get_sample_names(cachexia)
  cachexia = set_metadata(cachexia, new.metadata)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ metadata }
\keyword{ dataset }% __ONLY ONE__ keyword per line
