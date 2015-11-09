\name{pca_plot_3d}
\alias{pca_plot_3d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  3D pca plot
}
\description{
  3D plot from 3 components
}
\usage{
pca_plot_3d(dataset, model, var.class, pcas = 1:3, colors = NULL, 
legend.place = "topright", ...)
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
  data(cachexia)
  train.result = train_models_performance(cachexia, "pls", 
	"Muscle.loss", "cv")
  pca_plot_3d(cachexia, train.result$final.models$pls, "Muscle.loss")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ pca }
\keyword{ pls }% __ONLY ONE__ keyword per line
