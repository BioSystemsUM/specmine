\name{convert_from_chemospec}
\alias{convert_from_chemospec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Convert from ChemoSpec
}
\description{
   Convert the dataset in the ChemoSpec format to a dataset of this package.
}
\usage{
convert_from_chemospec(csobj, type = "undefined", 
description = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{csobj}{
   ChemoSpec object representing the dataset.
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
\dontrun{
  ## Example of converting a dataset from ChemoSpec 
  dataset = convert_from_chemospec(chemospec.dataset, type = "nmr-peaks", 
   	    description = "some description") 
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ chemospec }
\keyword{ dataset }% __ONLY ONE__ keyword per line
