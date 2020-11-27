\name{aov_all_vars}
\alias{aov_all_vars}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Analysis of variance
}
\description{
  Perform analysis of variance of all variables in the dataset.
}
\usage{
aov_all_vars(dataset, column.class, doTukey = TRUE, write.file = FALSE, 
file.out = "anova-res.csv")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{column.class}{
    string or index indicating what metadata to use.
}
  \item{doTukey}{
    boolean value for do or do not TukeyHSD.
}
  \item{write.file}{
    boolean value indicating if a file with the results is written or not.  
}
  \item{file.out}{
    name of the file if write.file is TRUE.
}
}

\value{
Data frame with the results of ANOVA, with p-value, logarithm of p-value, false discovery rate (fdr) and tukey is doTukey is TRUE. The result is ordered by p-value. 
}

\examples{
  ## Example of ANOVA with TukeyHSD
  library(specmine.datasets)
  data(propolis)
  propolis_proc = missingvalues_imputation(propolis)
  propolis_proc = flat_pattern_filter(propolis_proc, "iqr", by.percent = TRUE,
	red.value = 75) 
  result = aov_all_vars(propolis_proc, "seasons", doTukey = FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ anova }
\keyword{ tukey }% __ONLY ONE__ keyword per line
