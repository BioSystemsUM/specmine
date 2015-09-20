\name{find.equal.samples}
\alias{find.equal.samples}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Find equal samples
}
\description{
  Finds samples that have the same peak values - x and y (equal data frames)
}
\usage{
find.equal.samples(sample.list)
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
\dontrun{
    ## Example of finding equal samples
    equal.samples = find.equal.samples(sample.list)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~equal }
\keyword{ ~sample }% __ONLY ONE__ keyword per line
