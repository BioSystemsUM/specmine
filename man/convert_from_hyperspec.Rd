\name{convert_from_hyperspec}
\alias{convert_from_hyperspec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
   Convert from hyperspec
}
\description{
   Convert the dataset in the hyperspec format to a dataset of this package.
}
\usage{
convert_from_hyperspec(hsobj, type = "undefined", 
description = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hsobj}{
   hyperspec object representing the dataset.
}
  \item{type}{
   string representing the type of the data.
}
  \item{description}{
   string representing the description of the dataset.
}
}

\value{
   Returns a list representing the dataset converted.
}

\examples{
  ## Example of converting a dataset from hyperspec 
  library(hyperSpec)
  data(flu)
  dataset = convert_from_hyperspec(flu, type = "undefined", 
	    description = "some description") 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hyperspec }
\keyword{ dataset }% __ONLY ONE__ keyword per line
