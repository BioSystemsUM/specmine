\name{aov.one.var}
\alias{aov.one.var}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Analysis of variance of one variable
}
\description{
  Perform analysis of variance of one variable.
}
\usage{
aov.one.var(dataset, x.val, groups, doTukey = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{x.val}{
    variable name.
}
  \item{groups}{
    vector with samples groups.
}
  \item{doTukey}{
    boolean value for do or do not TukeyHSD.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
List with the results:
  \item{pvalue }{p-value of the variable}
  \item{log }{the logarithm of p-value}
  \item{fdr }{false discovery rate of the variable}
  \item{tukey }{results of TukeyHSD only if doTukey is TRUE}
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
  ## Example of anova one variable with TukeyHSD
  result = aov.one.var(metabolomics.dataset, "var.name", samples.groups, doTukey = TRUE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~anova }
\keyword{ ~tukey }% __ONLY ONE__ keyword per line
