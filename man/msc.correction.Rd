\name{msc.correction}
\alias{msc.correction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multiplicative scatter correction
}
\description{
  Perform multiplicative scatter correction on the spectra.
}
\usage{
msc.correction(dataset)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  Return the dataset with the multiplicative scatter correction employed on the data.
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
  ## Example of multiplicative scatter correction
  dataset.msc = msc.correction(dataset)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~msc }
\keyword{ ~dataset }% __ONLY ONE__ keyword per line
