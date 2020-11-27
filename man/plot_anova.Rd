\name{plot_anova}
\alias{plot_anova}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
   Plot ANOVA results
}
\description{
   Function for plotting the results from ANOVA.
}
\usage{
plot_anova(dataset, anova.results, anova.threshold = 0.01, 
reverse.x = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{anova.results}{
    ANOVA results.
}
  \item{anova.threshold}{
    ANOVA threshold for the p-value.
}
  \item{reverse.x}{
    boolean value to indicate if the x-axis is plotted in reverse.
}
}



\examples{
    ## Example of plotting the ANOVA results - first filter the 
    ## dataset to reduce computation time
	library(specmine.datasets)
	data(propolis)
    propolis_proc = missingvalues_imputation(propolis)
    propolis_proc = flat_pattern_filter(propolis_proc, "iqr", by.percent = TRUE,
	red.value = 75) 
    anova.results = aov_all_vars(propolis_proc, "seasons", doTukey = FALSE)
    plot_anova(propolis_proc, anova.results)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ anova }
\keyword{ plot }% __ONLY ONE__ keyword per line
