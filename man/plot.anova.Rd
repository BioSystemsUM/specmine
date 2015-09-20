\name{plot.anova}
\alias{plot.anova}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
   Plot ANOVA results
}
\description{
   Function for plotting the results from ANOVA.
}
\usage{
plot.anova(dataset, anova.results, anova.threshold = 0.01, reverse.x = F)
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
  \dontrun{
    ## Example of plotting the ANOVA results
    plot.anova(dataset, anova.results, 0.05)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~anova }
\keyword{ ~plot }% __ONLY ONE__ keyword per line
