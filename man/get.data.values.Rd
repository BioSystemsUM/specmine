\name{get.data.values}
\alias{get.data.values}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get data values
}
\description{
Gets the values of all samples in the dataset given a set of x axis names or indexes.
}
\usage{
get.data.values(dataset, x.axis.val, by.index = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
list representing the dataset from a metabolomics experiment.
}
  \item{x.axis.val}{
vector with the values of the x axis (could be names or indexes).
}
  \item{by.index}{
boolean value indicating if the x.axis.val is a vector of indexes or not.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a matrix with the values of all samples in the specified x axis.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
\dontrun{
  ## Example of getting a metadata value
  data.values = get.data.values(dataset, c("var1","var3","var4"), by.index = F)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~data }
\keyword{ ~values }% __ONLY ONE__ keyword per line
