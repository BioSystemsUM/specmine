\name{predict_samples}
\alias{predict_samples}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Predict samples
}
\description{
  Predict new samples.
}
\usage{
predict_samples(train.result, new.samples)
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

\value{
  Returns a data frame with the samples and the predicted class.
}

\examples{
  ## Example of predicting samples
  #data(cachexia)
  #training.result = train_models_performance(cachexia, "pls", 
  #	"Muscle.loss", "cv")
  #result = predict_samples(training.result$final.models$pls, cachexia$data)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ predict }
\keyword{ samples }% __ONLY ONE__ keyword per line
