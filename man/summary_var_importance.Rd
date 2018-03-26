\name{summary_var_importance}
\alias{summary_var_importance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Summary of variables importance
}
\description{
  Summary of variables importance of the models
}
\usage{
summary_var_importance(performances, number.rows)
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

\value{
  Returns list with the variables importance of each model.
}

\examples{
  ## Example of getting a summary of variables importance
  #data(cachexia)
  #training.result = train_models_performance(cachexia, "pls", 
  #	"Muscle.loss", "cv")
  #result = summary_var_importance(training.result, 10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ vip }
\keyword{ summary }% __ONLY ONE__ keyword per line
