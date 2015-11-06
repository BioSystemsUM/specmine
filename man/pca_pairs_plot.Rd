\name{pca_pairs_plot}
\alias{pca_pairs_plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  PCA pairs plot
}
\description{
  Shows a PCA pairs plot.
}
\usage{
pca_pairs_plot(dataset, pca.result, column.class = NULL, 
pcas = c(1, 2, 3, 4, 5), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
        list representing the dataset from a metabolomics experiment.
}
  \item{pca.result}{
    prcomp object with the PCA results.
}
  \item{column.class}{
    metadata's variable.
}
  \item{pcas}{
    the principal components to be shown.
}
  \item{\dots}{
	additional parameters to ggpairs function from GGally package.
}
}

\examples{
  ## Example of a PCA pairs plot
  data(cachexia)
  pca.result = pca_analysis_dataset(cachexia)
  pca_pairs_plot(cachexia, pca.result, "Muscle.loss", pcas = c(1,2,3))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ pairs }
\keyword{ pca }% __ONLY ONE__ keyword per line
