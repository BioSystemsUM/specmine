\name{find_equal_samples}
\alias{find_equal_samples}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Find equal samples
}
\description{
  Finds samples that have the same peak values - x and y (equal data frames)
}
\usage{
find_equal_samples(sample.list)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sample.list}{
        list of data frames with the samples' peaks.
}
}

\value{
  Returns a dataframe with two columns indicating which pair of samples are equal.
}

\examples{
    ## Example of finding equal samples
    data(propolisSampleList)
    equal.samples = find_equal_samples(propolisSampleList)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ equal }
\keyword{ sample }% __ONLY ONE__ keyword per line
