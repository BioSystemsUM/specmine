\name{train_and_predict}
\alias{train_and_predict}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Train and predict
}
\description{
  Train a model and predict new unlabeled samples with that model.
}
\usage{
train_and_predict(dataset, new.samples, column.class, model, 
validation, num.folds = 10, num.repeats = 10, tunelength = 10, 
tunegrid = NULL, metric = NULL, summary.function = 
defaultSummary)
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
  \item{metric}{
	metric used to evaluate the model's performance. Can be "Accuracy" or "ROC".
}
  \item{summary.function}{
	summary function. For "ROC" the multiClassSummary function must be used.
}
}

\value{
  Returns a list with the training result and the predictions result.
}

\examples{
  ## Example of training and predicting
  library(specmine.datasets)
  data(cachexia)
  result = train_and_predict(cachexia, new.samples = cachexia$data,
	 "Muscle.loss", "pls", "cv")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ train }
\keyword{ predict }% __ONLY ONE__ keyword per line
