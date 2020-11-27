\name{plot_fold_change}
\alias{plot_fold_change}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
   Plot fold change results
}
\description{
   Function for plotting the results from fold change.
}
\usage{
plot_fold_change(dataset, fc.results, fc.threshold, plot.log = TRUE, 
var = FALSE, xlab = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{fc.results}{
    fold change results.
}
  \item{fc.threshold}{
    fold change threshold for the p-value.
}
  \item{plot.log}{
    boolean value to determine if the fold change values are transformed logarithmically or not.
}
  \item{var}{
    boolean value, if TRUE it uses the xlab argument to represent the xlabel of the plot.
}
  \item{xlab}{
    string with the x axis description.
}
}

\examples{
    ## Example of plotting the fold change results
	library(specmine.datasets)
    data(cachexia)
    fc.results = fold_change(cachexia, "Muscle.loss", 
			"control")
    plot_fold_change(cachexia, fc.results, 2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ foldchange }
\keyword{ plot }% __ONLY ONE__ keyword per line
