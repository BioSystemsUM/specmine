\name{plotvar_twofactor}
\alias{plotvar_twofactor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Plot variable distribution on two factors
}
\description{
  Plot variable distribution on two factors from the dataset.
}
\usage{
plotvar_twofactor(dataset, variable, meta.var1, meta.var2, 
colour = "darkblue", title = "", xlabel = NULL, ylabel = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{variable}{
   variable's name.
}
  \item{meta.var1}{
   first metadata's variable.
}
  \item{meta.var2}{
   second metadata's variable.
}
  \item{colour}{
    colours of the bodies of the boxplots.
}
  \item{title}{
    title of the plot.
}
  \item{xlabel}{
   x-axis label.
}
  \item{ylabel}{
   y-axis label.
}
}

\value{
 Returns the plot object from ggplot function.
}

\examples{
  ## Example of plotting a variable's distribution with 2 factors
  data(cassavaPPD)
  plotvar_twofactor(cassavaPPD, "399.3", "varieties", "ppds")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ plot }
\keyword{ factor }% __ONLY ONE__ keyword per line
