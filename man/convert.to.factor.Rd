\name{convert.to.factor}
\alias{convert.to.factor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Convert metadata to factor
}
\description{
  Convert a metadata's variable to factor.
}
\usage{
convert.to.factor(dataset, metadata.var)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
  list representing the dataset from a metabolomics experiment.
}
  \item{metadata.var}{
  name of the metadata's variable.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  Returns the dataset with the metadata's variable converted to factor.
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
  ## Example of converting a metadata's variable to factor
  dataset = convert.to.factor(dataset, "metadata.var2")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~metadata }
\keyword{ ~factor }% __ONLY ONE__ keyword per line
