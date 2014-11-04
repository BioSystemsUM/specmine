\name{multifactor.aov.all.vars}
\alias{multifactor.aov.all.vars}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Multifactor ANOVA 
}
\description{
  Perform ANOVA on all variables with various metadata variables.
}
\usage{
multifactor.aov.all.vars(dataset, metadata.vars, write.file = F, file.out = "anova-res.csv")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{metadata.vars}{
    metadata variables to use in ANOVA.
}
  \item{write.file}{
    boolean value to write or not a file with the results
}
  \item{file.out}{
    name of the file.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  table with the p-values of all variables in the metadata variables.
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
  ## Example of multifactor ANOVA on all variables 
  pvalues = multifactor.aov.onevar(metabolomics.dataset, c("metadata1","metadata2"), write.file = T, file.out = "anova.csv")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~anova }
\keyword{ ~multifactor }% __ONLY ONE__ keyword per line
