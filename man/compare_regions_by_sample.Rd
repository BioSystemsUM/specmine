\name{compare_regions_by_sample}
\alias{compare_regions_by_sample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Compare regions by sample
}
\description{
 Compare two regions of a dataset by samples.
}
\usage{
compare_regions_by_sample(dataset1, dataset2, fn.to.apply, 
samples = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset1}{
    list representing the dataset from a metabolomics experiment.
}
  \item{dataset2}{
    list representing the dataset from a metabolomics experiment.
}
  \item{fn.to.apply}{
    function to apply (e.g. mean, max, min).
}
  \item{samples}{
    if defined restricts the application to a given set of samples.
}
  \item{\dots}{
    additional parameters to apply.by.sample function.
}
}

\value{
 Returns a data.frame with the results of the function applied to the samples and the ration between the two datasets.
}

\examples{
     ## Example of comparing regions by sample
	 library(specmine.datasets)
     data(cachexia)
     subset1 = subset_x_values(cachexia, 1:31, by.index = TRUE)
     subset2 = subset_x_values(cachexia, 32:63, by.index = TRUE)
     comp.regions.result = compare_regions_by_sample(subset1, subset2, 
			   mean)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ region }
\keyword{ sample }% __ONLY ONE__ keyword per line
