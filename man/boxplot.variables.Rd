\name{boxplot.variables}
\alias{boxplot.variables}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Boxplot of variables
}
\description{
  Boxplot of each variable of the dataset.
}
\usage{
boxplot.variables(dataset, variables = NULL, samples = NULL, horizontal = T, col = "lightblue", nchar.label = 10, cex.axis = 0.8, ...)
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
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
\dontrun{
  ## Example of showing the boxplot of a few variables
  boxplot.variables(dataset, variables = c("var1","var2","var3"))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~boxplot }
\keyword{ ~variables }% __ONLY ONE__ keyword per line
