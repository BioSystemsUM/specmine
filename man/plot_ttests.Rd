\name{plot_ttests}
\alias{plot_ttests}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
   Plot t-tests results
}
\description{
   Function for plotting the results from t-tests.
}
\usage{
plot_ttests(dataset, tt.results, tt.threshold = 0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{tt.results}{
    t-tests results.
}
  \item{tt.threshold}{
    t-test threshold for the p-value.
}
}



\examples{
    ## Example of plotting the t-tests results
    data(cachexia)
    ttests.results = tTests_dataset(cachexia, "Muscle.loss")
    plot_ttests(cachexia, ttests.results, 0.05)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ttest }
\keyword{ plot }% __ONLY ONE__ keyword per line
