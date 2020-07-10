\name{get_data_value}
\alias{get_data_value}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Get data value
}
\description{
  Get a data value given the x-axis labels and the sample
}
\usage{
get_data_value(dataset, x.axis.val, sample, by.index = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{x.axis.val}{
    index or name of the x-axis value.
}
  \item{sample}{
    index or name of the sample.
}
  \item{by.index}{
    boolean value indicating if the x-axis value and sample are represented as index or not.
}
}

\value{
  Returns a numeric with the data point value. 
}

\examples{
  ## Example of getting a data value from the dataset
  data(cachexia)
  data.value = get_data_value(cachexia, "Creatine", "PIF_178", 
	       by.index = FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ dataset }
\keyword{ value }% __ONLY ONE__ keyword per line
