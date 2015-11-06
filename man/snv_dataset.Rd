\name{snv_dataset}
\alias{snv_dataset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Standard Normal Variate
}
\description{
  Performs Standard Normal Variate on the dataset.
}
\usage{
snv_dataset(dataset)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
}

\value{
  Returns the dataset with the data normalized by SNV.
}


\examples{
  ## Example of SNV on a dataset
  data(cassavaPPD)
  dataset = snv_dataset(cassavaPPD)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ snv }
\keyword{ normalization }% __ONLY ONE__ keyword per line
