\name{pca.plot.3d}
\alias{pca.plot.3d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  3D pca plot
}
\description{
  3D plot from 3 components
}
\usage{
pca.plot.3d(dataset, model, var.class, pcas = 1:3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{model}{
    model with pca scores (pls model).
}
  \item{var.class}{
    metadata column class.
}
  \item{pcas}{
    the components to be plotted.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
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
  ## Example of a 3d pca plot
  pca.plot.3d(dataset, train.result$final.models$pls, "metadata_var")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~pca }
\keyword{ ~pls }% __ONLY ONE__ keyword per line
