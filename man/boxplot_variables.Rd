\name{boxplot_variables}
\alias{boxplot_variables}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Boxplot of variables
}
\description{
  Boxplot of each variable of the dataset.
}
\usage{
boxplot_variables(dataset, variables = NULL, samples = NULL, 
horizontal = TRUE, col = "lightblue", nchar.label = 10, 
cex.axis = 0.8, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{variables}{
    vector with the variables names or a NULL value indicating all variables.
}
  \item{samples}{
    vector with the samples names or a NULL value indicating all samples.
}
  \item{horizontal}{
    boolean value indicating if the boxplots should be horizontal.
}
  \item{col}{
    string that represents the color of the bodies of the boxplots.
}
  \item{nchar.label}{
    number of characters to display the variables' names.
}
  \item{cex.axis}{
    numeric value that indicates the amount by which the axis is magnified relative to the default.
}
  \item{\dots}{
  additional parameters of boxplot function.
}
}

\examples{
  ## Example of showing the boxplot of a few variables
  library(specmine.datasets)
  data(cachexia)
  boxplot_variables(cachexia, variables = c("Creatine","Serine","Lactate"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ boxplot }
\keyword{ variables }% __ONLY ONE__ keyword per line
