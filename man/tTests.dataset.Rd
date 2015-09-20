\name{tTests.dataset}
\alias{tTests.dataset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  t-Tests on dataset
}
\description{
  Run t-Tests for each row of the data from the dataset.
}
\usage{
tTests.dataset(dataset, metadata.var, threshold = NULL, write.file = F, file.out = "ttests.csv")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{metadata.var}{
    metadata variable to use in the t-tests.
}
  \item{threshold}{
    threshold value of the p-value.
}
  \item{write.file}{
    boolean value to write or not a file with the results.
}
  \item{file.out}{
    name of the file.
}
}

\value{
  Table with the results of the t-tests, with p-value, -log10(p-value) and false discovery rate (fdr).
}


\examples{
\dontrun{
  ## Example of t-Tests on dataset
  ttests.result = tTests.dataset(metabolomics.dataset, "class", write.file = T, file.out = "ttests.csv")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~ttest }
