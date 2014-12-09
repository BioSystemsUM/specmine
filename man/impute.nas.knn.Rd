\name{impute.nas.knn}
\alias{impute.nas.knn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Impute missing values with KNN
}
\description{
  Impute missing values with KNN
}
\usage{
impute.nas.knn(dataset, k = 10, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
      list representing the dataset from a metabolomics experiment.
}
  \item{k}{
    the number of nearest neighbors.
}
  \item{\dots}{
    additional values to impute.knn function.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
   Returns the dataset with no missing values.
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

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~missing }
\keyword{ ~values }% __ONLY ONE__ keyword per line
