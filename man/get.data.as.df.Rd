\name{get.data.as.df}
\alias{get.data.as.df}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Get data as data frame
}
\description{
  Get the data matrix from the dataset as a data frame.
}
\usage{
get.data.as.df(dataset)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  Returns the data matrix from the dataset as a data.frame object.
}

\examples{
\dontrun{
  ## Example of getting the data matrix as data frame
  datamatrix.dataframe = get.data.as.df(dataset)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~matrix }
\keyword{ ~dataframe }% __ONLY ONE__ keyword per line
