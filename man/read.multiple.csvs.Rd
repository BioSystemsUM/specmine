\name{read.multiple.csvs}
\alias{read.multiple.csvs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read multiple CSVs
}
\description{
  Reads multiple CSVs, each one with a sample.
}
\usage{
read.multiple.csvs(filenames, ext = ".csv", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filenames}{
    list of file names of the files to read.
}
  \item{ext}{
    extension name.
}
  \item{\dots}{
    additional parameters to read.csv function.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  returns a list of dataframes.
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
    filenames = c("sample1", "sample2", "sample3", "sample4")
    sample.list = read.multiple.csvs(filenames, ext = ".csv")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~read }
\keyword{ ~csv }% __ONLY ONE__ keyword per line
