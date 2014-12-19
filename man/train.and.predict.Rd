\name{train.and.predict}
\alias{train.and.predict}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Train and predict
}
\description{
  Train a model and predict new unlabeled samples with that model.
}
\usage{
train.and.predict(dataset, new.samples, column.class, model, validation, num.folds = 10, num.repeats = 10, tunelength = 10, tunegrid = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{new.samples}{
    dataframe with new samples to predict the class label.
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
  Returns a list with the training result and the predictions result.
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
  ## Example of training and predicting
  result = train.and.predict(dataset, new.samples, "metadata_var", "pls", "cv")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~train }
\keyword{ ~predict }% __ONLY ONE__ keyword per line
