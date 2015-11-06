\name{plot_kruskaltest}
\alias{plot_kruskaltest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
   Plot Kruskal-Wallis tests results
}
\description{
   Function for plotting the results from Kruskal-Wallis tests.
}
\usage{
plot_kruskaltest(dataset, kr.results, kr.threshold = 0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{kr.results}{
    Kruskal-Wallis tests results.
}
  \item{kr.threshold}{
    Kruskal-Wallis test threshold for the p-value.
}
}



\examples{
    ## Example of plotting the Kolmogorov-Smirnov tests results
    data(cachexia)
    kr.results = kruskalTest_dataset(cachexia, "Muscle.loss", 
	write.file = FALSE)
    plot_kruskaltest(cachexia, kr.results, 0.05)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Kruskal-Wallis }
\keyword{ plot }% __ONLY ONE__ keyword per line
