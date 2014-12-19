\name{train.classifier}
\alias{train.classifier}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
train.classifier(dataset, column.class, model, validation, num.folds = 10, num.repeats = 10, tunelength = 10, tunegrid = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
      list representing the dataset from a metabolomics experiment.
}
  \item{column.class}{
    metadata column class.
}
  \item{model}{
    model to be used in training.
}
  \item{validation}{
    validation method.
}
  \item{num.folds}{
    number of folds in cross validation.
}
  \item{num.repeats}{
    number of repeats.
}
  \item{tunelength}{
    number of levels for each tuning parameters.
}
  \item{tunegrid}{
    dataframe with possible tuning values.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  Returns the train result object from caret.
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
  ## Example of training a classifier
  result = train.classifier(dataset, "metadata_var", "pls", "cv")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~train }
\keyword{ ~classifier }% __ONLY ONE__ keyword per line
