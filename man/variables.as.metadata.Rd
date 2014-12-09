\name{variables.as.metadata}
\alias{variables.as.metadata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Variables as metadata
}
\description{
  Use one or more data variables as metadata variables.
}
\usage{
  variables.as.metadata(dataset, variables, by.index = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{variables}{
    name or index of the variables that are going to be used as metadata variables.
}
  \item{by.index}{
    boolean value indicating if the variables are indexes or names
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  Returns the dataset with the variables removed from the data and added on the metadata.
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
  ## Example of using a variable as metadata variable
  dataset = variables.as.metadata(dataset, "Var 1", by.index = F)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~metadata }
\keyword{ ~variable }% __ONLY ONE__ keyword per line
