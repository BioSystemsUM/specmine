\name{read.metadata}
\alias{read.metadata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Read metadata
}
\description{
  Read the metadata from a file.
}
\usage{
read.metadata(filename, header.col = T, header.row = T, sep = ",")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{
  name of the file with the metadata.
}
  \item{header.col}{
    boolean value indicating if the metadata CSV file contains a header column with the name of the metadata variables.
}
  \item{header.row}{
  boolean value indicating if the metadata CSV file contains a header row with the name of the samples.
}
  \item{sep}{
  the separator character.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  Returns a data frame with the metadata.
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
  ## Example of reading metadata from a file
  metadata = read.metadata("metadata.csv")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~metadata }
\keyword{ ~csv }% __ONLY ONE__ keyword per line
