\name{pca_analysis_dataset}
\alias{pca_analysis_dataset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  PCA analysis (classical)
}
\description{
  Performs a classical PCA analysis over the dataset.
}
\usage{
pca_analysis_dataset(dataset, scale = TRUE, center = TRUE, 
write.file = FALSE, file.out = "pca", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{scale}{
    boolean value indicating if the variables are going to be scaled or not.
}
  \item{center}{
    booleam value indicating if the variables are going to be centered or not.
}
  \item{write.file}{
    boolean value that indicates if the results from PCA analysis are going to be written on a file.
}
  \item{file.out}{
    name of the file that will store the results.
}
  \item{\dots}{
    additional parameters to ggplot function.
}
}

\value{
  object of class 'prcomp' with the results from the PCA analysis.
}

\examples{
  ## Example of performing a classical PCA analysis
  data(cachexia)
  pca.results = pca_analysis_dataset(cachexia)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ pca }
\keyword{ unsupervised }% __ONLY ONE__ keyword per line
