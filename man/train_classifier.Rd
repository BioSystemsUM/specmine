\name{train_classifier}
\alias{train_classifier}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Train classifier
}
\description{
  Train a specific classifier.
}
\usage{
train_classifier(dataset, column.class, model, validation, 
num.folds = 10, num.repeats = 10, tunelength = 10, 
tunegrid = NULL, metric = NULL, 
summary.function = defaultSummary, class.in.metadata = T)
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
  \item{metric}{
	metric used to evaluate the model's performance. Can be "Accuracy" or "ROC".
}
  \item{summary.function}{
	summary function. For "ROC" the multiClassSummary function must be used.
}
  \item{class.in.metadata}{
	boolean value to indicate if the class is in metadata.
}
}

\value{
  Returns the train result object from caret.
}

\examples{
  ## Example of training a classifier
  data(cachexia)
  result = train_classifier(cachexia, "Muscle.loss", "pls", "cv")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ train }
\keyword{ classifier }% __ONLY ONE__ keyword per line
