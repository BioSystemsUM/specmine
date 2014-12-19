\name{predict.samples}
\alias{predict.samples}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Predict samples
}
\description{
  Predict new samples.
}
\usage{
predict.samples(train.result, new.samples)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{train.result}{
    result from training a classifier.
}
  \item{new.samples}{
    dataframe with new samples.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  Returns a data frame with the samples and the predicted class.
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
  ## Example of predicting samples
  result = predict.samples(training.result, new.samples)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~predict }
\keyword{ ~samples }% __ONLY ONE__ keyword per line
