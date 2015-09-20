\name{pca.plot.3d}
\alias{pca.plot.3d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  3D pca plot
}
\description{
  3D plot from 3 components
}
\usage{
pca.plot.3d(dataset, model, var.class, pcas = 1:3, colors = NULL, legend.place = "topright", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{model}{
    model with pca scores (pls model).
}
  \item{var.class}{
    metadata column class.
}
  \item{pcas}{
    the components to be plotted.
}
  \item{colors}{
	colors of the groups.
}
  \item{legend.place}{
	legend place.
}
  \item{\dots}{
	additional parameters to legend function.
}
}

\examples{
\dontrun{
  ## Example of a 3d pca plot
  pca.plot.3d(dataset, train.result$final.models$pls, "metadata_var")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~pca }
\keyword{ ~pls }% __ONLY ONE__ keyword per line
