\name{train_models_performance}
\alias{train_models_performance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Train models
}
\description{
  Train various models.
}
\usage{
train_models_performance(dataset, models, column.class, 
validation, num.folds = 10, num.repeats = 10, tunelength = 10, 
tunegrid = NULL, metric = NULL, summary.function = "default", 
class.in.metadata = T, compute.varimp = T)
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
  \item{metric}{
	metric used to evaluate the model's performance. Can be "Accuracy" or "ROC".
}
  \item{summary.function}{
	summary function. For "ROC" the multiClassSummary function must be used.
}
  \item{class.in.metadata}{
	boolean value to indicate if the class is in metadata.
}
  \item{compute.varimp}{
	boolean value to indicate if the var importance is calculated.
}
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

\examples{
  ## Example of training models
  #data(cachexia)
  #result = train_models_performance(cachexia, "pls", 
  #	   "Muscle.loss", "cv")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ train}
\keyword{ model }% __ONLY ONE__ keyword per line
