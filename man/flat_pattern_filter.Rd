\name{flat_pattern_filter}
\alias{flat_pattern_filter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Flat pattern filter
}
\description{
  Performs a flat pattern filter over the dataset.
}
\usage{
flat_pattern_filter(dataset, filter.function = "iqr", 
by.percent = TRUE, by.threshold = FALSE, red.value = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{filter.function}{
    filter function. It can be:
    \itemize{
      \item \bold{"iqr"} - Interquantile Range
      \item \bold{"rsd"} - Relative Standard Deviation
      \item \bold{"sd"} - Standard Deviation
      \item \bold{"mad"} - Median Absolute Deviation
      \item \bold{"mean"} - Mean
      \item \bold{"median"} - Median
    }
}
  \item{by.percent}{
    boolean value, if TRUE the number of variables to filter will be a percentage of the number of variables in the dataset; percentage is given by the "red.value" parameter
}
  \item{by.threshold}{
    boolean value, if TRUE, defines filtering will select variables where values of filtering function are below a given threshold. Threshold is defined by red.value that defines the minimum value of the function needed to keep the variable.
}
  \item{red.value}{
    it can be the percentage or the threshold number. If red.value = "auto", will calculate number of variables to remove automatically
}
}

\value{
  Returns the dataset with the data filtered.
}
\examples{
  ## Example of flat pattern filter
  library(specmine.datasets)
  data(propolis)
  dataset.filtered = flat_pattern_filter(propolis, "iqr", by.percent = TRUE, 
		     red.value = 20)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ flatpattern }
\keyword{ filter }% __ONLY ONE__ keyword per line
