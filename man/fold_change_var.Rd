\name{fold_change_var}
\alias{fold_change_var}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Fold change applied on two variables
}
\description{
 Fold change applied on two variables. Instead of having the difference of the variables on two groups, we have the difference of the groups on two variables.
}
\usage{
fold_change_var(dataset, metadata.var, variables, 
threshold.min.fc = NULL, write.file = FALSE, 
file.out = "fold_change_reverse.csv")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{metadata.var}{
    metadata to use to calculate the fold change.
}
  \item{variables}{
    vector with two positions containing the name of the variables.
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
  ## Example of fold change reverse
  library(specmine.datasets)
  data(cachexia)
  fold.change.results = fold_change_var(cachexia, "Muscle.loss", 
			c("Creatine", "Serine"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ foldchange }
\keyword{ reverse }% __ONLY ONE__ keyword per line
