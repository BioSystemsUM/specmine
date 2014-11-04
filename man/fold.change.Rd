\name{fold.change}
\alias{fold.change}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Fold change analysis
}
\description{
  Perform fold change analysis on the dataset.
}
\usage{
fold.change(dataset, metadata.var, ref.value, write.file = F, file.out = "fold_change.csv")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{metadata.var}{
    metadata to use to calculate the fold change.
}
  \item{ref.value}{
    class name to indicate the initial value.
}
  \item{write.file}{
    boolean value to write or not a file with the results.
}
  \item{file.out}{
    name of the file.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  table of results with fold change and log2 of fold change. 
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
  fold.change.results = fold.change(metabolomics.dataset, "metadata", "initial.value", write.file = T, file.out = "foldchange.csv")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~foldchange }
