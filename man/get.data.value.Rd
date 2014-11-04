\name{get.data.value}
\alias{get.data.value}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Get data value
}
\description{
  Get a data value given the x-axis labels and the sample
}
\usage{
get.data.value(dataset, x.axis.val, sample, by.index = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{x.axis.val}{
    index or name of the x-axis value.
}
  \item{sample}{
    index or name of the sample.
}
  \item{by.index}{
    boolean value indicating if the x-axis value and sample are represented as index or not.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  Returns a numeric with the data point value. 
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
\dontrun{
  ## Example of getting a data value from the dataset
  data.value = get.data.value(dataset, "xaxis.var1", "sample2", by.index = F)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~dataset }
\keyword{ ~value }% __ONLY ONE__ keyword per line
