\name{remove_samples_by_na_metadata}
\alias{remove_samples_by_na_metadata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Remove samples by NA on metadata
}
\description{
 Remove samples from the dataset with the metadata's variable value with NAs.
}
\usage{
remove_samples_by_na_metadata(dataset, metadata.var)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{metadata.var}{
    metadata's variable.
}
}

\value{
 Returns the dataset with the samples with NA on metadata's variable removed.
}

\examples{
  ## Example of removing samples that have NAs on metadata's variable
  library(specmine.datasets)
  data(cachexia)
  cachexia$metadata$Muscle.loss[10] = NA
  cachexia = remove_samples_by_na_metadata(cachexia, "Muscle.loss")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ remove }
\keyword{ metadata }% __ONLY ONE__ keyword per line
