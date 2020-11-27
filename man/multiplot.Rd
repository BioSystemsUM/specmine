\name{multiplot}
\alias{multiplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Multiplot
}
\description{
 Multiplot from ggplot2 - function taken from (see references).
}
\usage{
multiplot(plots, plotlist = NULL, file, cols = 1, layout = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{plots}{
   list with the plots to display.
}
  \item{plotlist}{
   plot list.
}
  \item{file}{
   file.
}
  \item{cols}{
   number of columns.
}
  \item{layout}{
   layout of the plot.
}
}

\references{
 http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_%28ggplot2%29/
}

\examples{
     ## Example of multiplot
	 library(specmine.datasets)
     data(cachexia)
     pca.result = pca_analysis_dataset(cachexia)
     plot1 = pca_scoresplot2D(cachexia, pca.result, "Muscle.loss", 
	ellipses = TRUE)
     plot2 = pca_scoresplot2D(cachexia, pca.result, "Muscle.loss", 
	ellipses = FALSE, labels = TRUE)
     plts = list(plot1,plot2)
     multiplot(plts, cols = 2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multiplot }
\keyword{ ggplot }% __ONLY ONE__ keyword per line
