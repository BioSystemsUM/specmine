\name{low_level_fusion}
\alias{low_level_fusion}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Low level fusion
}
\description{
Low level fusion method for integrate different datasets (only samples with the same name on all datasets will be merged)
}
\usage{
low_level_fusion(datasets)
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

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ fusion }
\keyword{ dataset }% __ONLY ONE__ keyword per line
