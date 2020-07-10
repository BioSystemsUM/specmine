\name{plot_regression_coefs_pvalues}
\alias{plot_regression_coefs_pvalues}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Plot regression coefficient and p-values
}
\description{
 Plots the linear regression coefficient and the p-values.
}
\usage{
plot_regression_coefs_pvalues(linreg.results, bar.col = NULL, 
coef.size = 5, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{linreg.results}{
   linear regression results.
}
  \item{bar.col}{
   color of the bars.
}
  \item{coef.size}{
   coefficient font size.
}
  \item{\dots}{
   additional parameters to geom_text and geom_bar from ggplot.
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
\keyword{ plot }% __ONLY ONE__ keyword per line
