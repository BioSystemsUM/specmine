\name{plot_kstest}
\alias{plot_kstest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
   Plot Kolmogorov-Smirnov tests results
}
\description{
   Function for plotting the results from Kolmogorov-Smirnov tests.
}
\usage{
plot_kstest(dataset, ks.results, ks.threshold = 0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{ks.results}{
    Kolmogorov-Smirnov tests results.
}
  \item{ks.threshold}{
    Kolmogorov-Smirnov test threshold for the p-value.
}
}



\examples{
    ## Example of plotting the Kolmogorov-Smirnov tests results
	library(specmine.datasets)
    data(cachexia)
    ks.results = ksTest_dataset(cachexia, "Muscle.loss", 
	write.file = FALSE)
    plot_kstest(cachexia, ks.results, 0.05)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Kolmogorov-Smirnov }
\keyword{ plot }% __ONLY ONE__ keyword per line
