\name{summary.var.importance}
\alias{summary.var.importance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Summary of variables importance
}
\description{
  Summary of variables importance of the models
}
\usage{
summary.var.importance(performances, number.rows)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{performances}{
    the result from training the models.
}
  \item{number.rows}{
    number of variables to display.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  Returns list with the variables importance of each model.
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
  ## Example of getting a summary of variables importance
  result = summary.var.importance(training.result, 10)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~vip }
\keyword{ ~summary }% __ONLY ONE__ keyword per line
