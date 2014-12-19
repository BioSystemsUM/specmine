\name{train.models.performance}
\alias{train.models.performance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Train models
}
\description{
  Train various models.
}
\usage{
train.models.performance(dataset, models, column.class, validation, num.folds = 10, num.repeats = 10, tunelength = 10, tunegrid = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
      list representing the dataset from a metabolomics experiment.
}
  \item{models}{
    models to be used in training.
}
  \item{column.class}{
    metadata column class.
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
Returns a list with the results from training
\item{performance}{The results from the best tunes of the models}
\item{vips}{The variable importance from the models}
\item{full.results}{The full results from the tuning parameters of each model}
\item{best.tunes}{The best tune of each model}
\item{confusion.matrices}{The confusion matrices of the models (only in classification)}
\item{final.models}{The final models}
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
  ## Example of training models
  result = train.models.performance(dataset, c("pls","rf"), "metadata_var", "cv")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~train}
\keyword{ ~model }% __ONLY ONE__ keyword per line
