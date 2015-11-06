\name{volcano_plot_fc_tt}
\alias{volcano_plot_fc_tt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
   Volcano plot
}
\description{
   Volcano plot to intersect the results from t-tests and fold change.
}
\usage{
volcano_plot_fc_tt(dataset, fc.results, tt.results, 
fc.threshold = 2, tt.threshold = 0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{fc.results}{
    fold change results.
}
  \item{tt.results}{
    t-tests results.
}
  \item{fc.threshold}{
    fold change threshold value.
}
  \item{tt.threshold}{
    t-test p-value threshold.
}
}

\value{
    Returns the name of the samples which intersects both fold change and t-tests results above the given thresholds.
}

\examples{
    ## Example of a volcano plot
    data(cachexia)
    foldchange.results = fold_change(cachexia, "Muscle.loss", "control")
    ttests.results = tTests_dataset(cachexia, "Muscle.loss")
    volcano_plot_fc_tt(cachexia, foldchange.results, ttests.results, 
		       fc.threshold = 2, tt.threshold = 0.01)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ volcano }
\keyword{ plot }% __ONLY ONE__ keyword per line
