\name{tTests.pvalue}
\alias{tTests.pvalue}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  t-Tests on matrix
}
\description{
  Run t-Tests on each row of the matrix.
}
\usage{
tTests.pvalue(datamat, samp.classes)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datamat}{
    matrix with the data.
}
  \item{samp.classes}{
    metadata variable to use in the t-tests.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  Data frame with p-values sorted in descending order.
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
  ## Example of t-Tests
  ttests.result = tTests.pvalue(matrix, "class")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~ttest }
