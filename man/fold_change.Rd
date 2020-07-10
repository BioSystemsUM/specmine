\name{fold_change}
\alias{fold_change}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Fold change analysis
}
\description{
  Perform fold change analysis on the dataset.
}
\usage{
fold_change(dataset, metadata.var, ref.value, 
threshold.min.fc = NULL, write.file = FALSE, 
file.out = "fold_change.csv")
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
  \item{threshold.min.fc}{
    minimum threshold of the fold change value.
}
  \item{write.file}{
    boolean value to write or not a file with the results.
}
  \item{file.out}{
    name of the file.
}
}

\value{
  Table of results with fold change and log2 of fold change. 
}


\examples{
  ## Example of fold change
  data(cachexia)
  fold.change.results = fold_change(cachexia, "Muscle.loss", 
			"control", write.file = FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ foldchange }
