\name{count.missing.values.per.sample}
\alias{count.missing.values.per.sample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Count missing values per sample
}
\description{
  Counts the missing values on each sample of the dataset.
}
\usage{
count.missing.values.per.sample(dataset, remove.zero = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{remove.zero}{
    boolean value indicating if the results of samples with no missing value are removed.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  Returns a vector with the number of missing values on each sample.
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
  ## Example of counting the missing values on each sample
  count.missing.values.per.sample(dataset)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~missing }
\keyword{ ~values }% __ONLY ONE__ keyword per line
