\name{xvalue.interval.to.indexes}
\alias{xvalue.interval.to.indexes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Get indexes of an interval of x-values
}
\description{
Returns indexes corresponding to an interval of x-values (only to numerical values - spectra)
}
\usage{
xvalue.interval.to.indexes(dataset, min.value, max.value)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
  list representing the dataset from a metabolomics experiment.
}
  \item{min.value}{
mininum x-value of the interval.
}
  \item{max.value}{
maximum x-value of the interval.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a numeric vector with the indexes of the x-values interval
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
  ## Example of getting a the indexes of an interval of x-values
  indexes.interval = xvalue.interval.to.indexes(dataset, 2.0, 5.5)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~indexes }
\keyword{ ~interval }% __ONLY ONE__ keyword per line
