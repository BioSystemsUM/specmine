\name{data_correction}
\alias{data_correction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Data correction
}
\description{
  Perform spectra corrections with 3 different methods.
}
\usage{
data_correction(dataset, type = "background", 
method = "modpolyfit", ...)
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

\value{
  Returns the dataset with the spectra corrected.
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ correction }
\keyword{ dataset }% __ONLY ONE__ keyword per line
