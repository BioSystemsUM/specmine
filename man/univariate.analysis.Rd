\name{univariate.analysis}
\alias{univariate.analysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Univariate Analysis
}
\description{
  Perform different methods of univariate analysis.
}
\usage{
univariate.analysis(dataset, type = "anova", column.class, ref.class, method = "pearson")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{type}{
    a string describing the type of univariate analysis. Possible values are:
    \itemize{
      \item "anova" perform analysis of variance.
      \item "ttest" perform t-tests.
      \item "foldchange" calculate fold-change between two classes.
      \item "correlation" calculate the correlations between the variables.
    }
}
  \item{column.class}{
    string or index indicating what metadata to use.
}
  \item{ref.class}{
    string indicating which class value is the initial value for fold change.
}
  \item{method}{
    string indicating the correlation method. Possible values are "pearson", "kendall" and "spearman".
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  List with the corresponding results of the method chosen.
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
  ## Example of fold change
  foldchange.result = univariate.analysis(metabolomics.dataset, type = "foldchange", "class.column", "initial.class")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~univariate }
\keyword{ ~anova }% __ONLY ONE__ keyword per line
\keyword{ ~foldchange }
\keyword{ ~correlations }
\keyword{ ~ttest }