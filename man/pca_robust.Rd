\name{pca_robust}
\alias{pca_robust}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  PCA analysis (robust)
}
\description{
  Performs a robust PCA analysis.
}
\usage{
pca_robust(dataset, center = "median", scale = "mad", k = 10, 
write.file = F, file.out = "robpca", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{center}{
    indicates how the data is to be centered. Can be a function or a vector with the center values of each column.
}
  \item{scale}{
 indicates how the data is to be rescaled. Can be a function or a vector with the scale value of each column.
}
  \item{k}{
    the desired number of components to compute
}
  \item{write.file}{
    boolean value that indicates if the results from PCA analysis are going to be written on a file.
}
  \item{file.out}{
    name of the file that will store the results.
}
  \item{\dots}{
    additional parameters pass to or from other functions.
}
}

\value{
  Returns an object of class 'princomp' with the PCA results.
}

\examples{
  ## Example of performing a robust PCA analysis
  data(cachexia)
  pca.results = pca_robust(cachexia, center = "mean", scale = "mad", 
		k = 10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ pca }
\keyword{ robust }% __ONLY ONE__ keyword per line
