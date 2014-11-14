\name{set.value.label}
\alias{set.value.label}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Set value label
}
\description{
  Set a new value label for the dataset.
}
\usage{
set.value.label(dataset, new.val.label)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
  list representing the dataset from a metabolomics experiment.
}
  \item{new.val.label}{
  string with the new value label.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  Returns the dataset with the updated value label.
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
  ## Example of setting a new value label to the dataset
  dataset = set.value.label(dataset, "new value label")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~label }
\keyword{ ~dataset }% __ONLY ONE__ keyword per line
