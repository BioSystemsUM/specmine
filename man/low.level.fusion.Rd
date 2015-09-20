\name{low.level.fusion}
\alias{low.level.fusion}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Low level fusion
}
\description{
Low level fusion method for integrate different datasets (only samples with the same name on all datasets will be merged)
}
\usage{
low.level.fusion(datasets)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datasets}{
 list containing the datasets to be fused (each dataset is a list from the pre defined format for the dataset).
}
}

\value{
Return a single dataset with all the datasets merged.
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
\dontrun{
  ## Example of low level fusion
  datasets = list(dataset1, dataset2)
  dataset.fusion = low.level.fusion(datasets)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~fusion }
\keyword{ ~dataset }% __ONLY ONE__ keyword per line
