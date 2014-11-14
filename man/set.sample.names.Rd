\name{set.sample.names}
\alias{set.sample.names}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Set samples names
}
\description{
  Set new samples names to the dataset.
}
\usage{
set.sample.names(dataset, new.sample.names)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
  list representing the dataset from a metabolomics experiment.
}
  \item{new.sample.names}{
  vector with the new samples names.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  Returns the dataset with the updated samples names.
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
  ## Example of setting a new value label to the dataset
  dataset = set.sample.names(dataset, new.samples.names)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~sample }
\keyword{ ~dataset }% __ONLY ONE__ keyword per line
