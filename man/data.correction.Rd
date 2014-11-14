\name{data.correction}
\alias{data.correction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Data correction
}
\description{
  Perform spectra corrections with 3 different methods.
}
\usage{
data.correction(dataset, type = "background", method = "modpolyfit", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{type}{
    string that represents the type of correction that will be applied to the spectra. The three possible types are: "background", to perform background correction; "offset", to perform offset correction; and "baseline", to perform baseline correction.
}
  \item{method}{
    string parameter of baseline correction indicating the correction method.
}
  \item{\dots}{
  additional parameters of baseline correction.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  Returns the dataset with the spectra corrected.
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
  ## Example of data correction using background correction
  dataset.corrected = data.correction(dataset, type = "background")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~correction }
\keyword{ ~dataset }% __ONLY ONE__ keyword per line
