\name{read.csvs.folder}
\alias{read.csvs.folder}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Read CSVs from folder
}
\description{
  Reads multiple CSV files in a given folder.
}
\usage{
read.csvs.folder(foldername, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{foldername}{
    string with the name of the folder.
}
  \item{\dots}{
    additional parameters to read.csv function.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  Returns a list of data frames.
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
    ## Example of reading multiple csvs
    sample.list = read.csvs.folder("foldername")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~read }
\keyword{ ~csv }% __ONLY ONE__ keyword per line
