\name{get.metadata.value}
\alias{get.metadata.value}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Get metadata value
}
\description{
  Get the metadata value
}
\usage{
get.metadata.value(dataset, variable, sample)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
list representing the dataset from a metabolomics experiment.
}
  \item{variable}{
index or name of the metadata variable.
}
  \item{sample}{
index or name of the sample.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Return the corresponding metadata value of the sample.
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
  metadata.value = get.metadata.value(dataset, "var1", "sample2")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
