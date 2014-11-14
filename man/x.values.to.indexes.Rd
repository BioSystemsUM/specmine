\name{x.values.to.indexes}
\alias{x.values.to.indexes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Get x-values indexes
}
\description{
  Returns the indexes corresponding to a vector of x-values (only to numerical values - spectra)
}
\usage{
x.values.to.indexes(dataset, x.values)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
list representing the dataset from a metabolomics experiment.
}
  \item{x.values}{
  vector of x-values.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  Returns a numeric vector with the indexes of the x-values.
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
  ## Example of getting the indexes of the x-values
  indexes = x.values.to.indexes(dataset, c(1.3, 3.5, 6))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~indexes }
\keyword{ ~xvalues }% __ONLY ONE__ keyword per line
