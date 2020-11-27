\name{heatmap_correlations}
\alias{heatmap_correlations}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Correlations heatmap
}
\description{
  Plots a heatmap with the correlations.
}
\usage{
heatmap_correlations(correlations, col = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{correlations}{
    correlation matrix
}
  \item{col}{
    colors to be used on heatmap.
}
  \item{\dots}{
    extra parameters to visual purposes.
}
}



\examples{
  ## Example of correlations heatmap
  library(specmine.datasets)
  data(cachexia)
  correlations = correlations_dataset(cachexia)
  heatmap_correlations(correlations)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ correlations }
\keyword{ heatmap }% __ONLY ONE__ keyword per line
