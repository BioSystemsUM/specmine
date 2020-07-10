\name{boxplot_vars_factor}
\alias{boxplot_vars_factor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Boxplot of variables with metadata's variable factors
}
\description{
 Boxplot of variables with metadata's variable factors from the dataset.
}
\usage{
boxplot_vars_factor(dataset, meta.var, variables = NULL, 
samples = NULL, horizontal = FALSE, nchar.label = 10, col = NULL, 
vec.par = NULL, cex.axis = 0.8, ylabs = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{meta.var}{
   metadata's variable.
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
  \item{nchar.label}{
    number of characters to display the variables' names.
}
  \item{col}{
    string that represents the color of the bodies of the boxplots.
}
  \item{vec.par}{
    vector with the disposition of the boxplots (rows, columns).
}
  \item{cex.axis}{
    numeric value that indicates the amount by which the axis is magnified relative to the default.
}
  \item{ylabs}{
    y-axis labels.
}
  \item{\dots}{
  additional parameters of boxplot function.
}
}


\examples{
  ## Example of showing the boxplot factors of a few variables
  data(cachexia)
  boxplot_vars_factor(cachexia, "Muscle.loss", variables = c("Creatine","Serine",
  "Lactate"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~boxplot }
\keyword{ ~factor }% __ONLY ONE__ keyword per line
