\name{mean_centering}
\alias{mean_centering}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Mean centering
}
\description{
 Performs mean centering on the dataset.
}
\usage{
mean_centering(dataset)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
     list representing the dataset from a metabolomics experiment.
}
}

\value{
 Returns the dataset with mean centering applied.
}

\examples{
  ## Example of mean centering
  data(cassavaPPD)
  cassavaPPD = mean_centering(cassavaPPD)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ mean }
\keyword{ centering }% __ONLY ONE__ keyword per line
