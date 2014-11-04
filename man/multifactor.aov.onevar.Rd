\name{multifactor.aov.onevar}
\alias{multifactor.aov.onevar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  One Variable Multifactor ANOVA
}
\description{
  Perform ANOVA on one variable with various metadata variables.
}
\usage{
multifactor.aov.onevar(dataset, x.val, metadata.vars)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{x.val}{
    variable to perform the ANOVA.
}
  \item{metadata.vars}{
    metadata variables to use in ANOVA.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  Returns the p-values of the variable in the metadata variables
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
  ## Example of multifactor ANOVA on one variable (1.1 ppm)
  pvalues = multifactor.aov.onevar(metabolomics.dataset, 1.1, c("metadata1","metadata2"))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~anova }
\keyword{ ~multifactor }% __ONLY ONE__ keyword per line
