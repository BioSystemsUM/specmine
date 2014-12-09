\name{metadata.as.variables}
\alias{metadata.as.variables}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Metadata as variables
}
\description{
  Use one or more metadata variables as variables.
}
\usage{
  metadata.as.variables(dataset, metadata.vars, by.index = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{metadata.vars}{
    name or index of the metadata variables that are going to be used as data  variables.
}
  \item{by.index}{
    boolean value indicating if the metadata variables are indexes or names
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  Returns the dataset with the metadata variables removed from the metadata and added on the data matrix.
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
  ## Example of using a metadata variable as data variable
  dataset = metadata.as.variables(dataset, "MetaVar 1", by.index = F)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~metadata }
\keyword{ ~variable }% __ONLY ONE__ keyword per line
