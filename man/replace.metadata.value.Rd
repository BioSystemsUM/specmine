\name{replace.metadata.value}
\alias{replace.metadata.value}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Replace metadata's value
}
\description{
  Replace a metadata's variable value of a sample.
}
\usage{
replace.metadata.value(dataset, variable, sample, new.value)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
  list representing the dataset from a metabolomics experiment.
}
  \item{variable}{
  metadata's variable.
}
  \item{sample}{
  name of the sample.
}
  \item{new.value}{
  new value of the metadata's variable.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  Returns the dataset with the metadata updated.
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
  ## Example of replacing metadata's variable value of a sample
  dataset = replace.metadata.value(dataset, "metadata.var1", "sample4", "new.value")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~metadata }
\keyword{ ~dataset }% __ONLY ONE__ keyword per line
