\name{set.x.values}
\alias{set.x.values}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Set new x-values
}
\description{
Set new x-values to the dataset
}
\usage{
set.x.values(dataset, new.x.values, new.x.label = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
  list representing the dataset from a metabolomics experiment.
}
  \item{new.x.values}{
  vector with the new x-values.
}
  \item{new.x.label}{
  string with the new x-label (can be NULL).
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  Returns the dataset with the updated x-values.
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
  ## Example of setting new x-values to the dataset
  dataset = set.x.values(dataset, new.xvalues, new.x.label = NULL)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~xvalues }
\keyword{ ~dataset }% __ONLY ONE__ keyword per line
