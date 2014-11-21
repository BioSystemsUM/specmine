\name{get.all.intensities}
\alias{get.all.intensities}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Get all intensities.
}
\description{
  Gets all intensities for a sample list; result is a data frame.
}
\usage{
get.all.intensities(sample.list, tol = 0.001)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sample.list}{
    list of data frames with the samples' peaks.
}
  \item{tol}{
    tolerance value of the frequencies.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  Returns a data frame with all the intensities of each sample in all frequencies.
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
    ## Example of getting the intensity
    df.intensities = get.all.intensities(sample, tol = 0.001)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~intensity }
\keyword{ ~frequencies }% __ONLY ONE__ keyword per line
