\name{calculate.shifts}
\alias{calculate.shifts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Calculate shifts
}
\description{
  Calculate shifts based on a band of the spectra.
}
\usage{
calculate.shifts(dataset, ref.limits = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{ref.limits}{
    vector with 2 elements that represents the reference limits to calculate the shifts.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  Returns a numeric vector with the shifts for each spectra.
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
  ## Example of calculating shifts
  shifts = calculate.shifts(dataset, ref.limits = c(200,250))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~shift }
\keyword{ ~spectra }% __ONLY ONE__ keyword per line
