\name{get_sample_2d_data}
\alias{get_sample_2d_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Get data
}
\description{
  Get the data matrix of a specific sample from a 2D dataset.
}
\usage{
get_sample_2d_data(dataset_2d, sample)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset_2d}{
  List representing the 2D dataset from a 2D metabolomics experiment.
  }
  
  \item{sample}{
  A string or integer, representing a sample from a 2D metabolomics experiment.
  }
}

\value{
  Returns the sample's data matrix.
}



% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ matrix }
\keyword{ 2d-dataset }% __ONLY ONE__ keyword per line