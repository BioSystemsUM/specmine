\name{set.metadata}
\alias{set.metadata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Set new metadata
}
\description{
Updates the dataset's metadata with a new one.
}
\usage{
set.metadata(dataset, new.metadata)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
  list representing the dataset from a metabolomics experiment.
}
  \item{new.metadata}{
  matrix or dataframe with the new metadata.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  Returns the dataset with the updated metadata.
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
  ## Example of setting a new metadata to the dataset
  dataset = set.metadata(dataset, new.metadata)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~metadata }
\keyword{ ~dataset }% __ONLY ONE__ keyword per line
