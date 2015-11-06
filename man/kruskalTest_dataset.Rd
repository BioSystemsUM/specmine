\name{kruskalTest_dataset}
\alias{kruskalTest_dataset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Kruskal-Wallis tests on dataset
}
\description{
  Run Kruskal-Wallis Tests for each row of the data from the dataset.
}
\usage{
kruskalTest_dataset(dataset, metadata.var, threshold = NULL, 
write.file = F, file.out = "kruskal.csv")
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
  Table with the results of the Kruskal-Wallis tests, with p-value, -log10(p-value) and false discovery rate (fdr).
}


\examples{
  ## Example of ks-Tests on dataset
  data(cachexia)
  kruskaltests.result = kruskalTest_dataset(cachexia, "Muscle.loss", 
		  write.file = FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Kruskal-Wallis }
