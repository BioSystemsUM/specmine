\name{find_corr}
\alias{find_corr}
\title{Find the best correlation value for the clustering of peaks}
\description{akes a dataset and calculates the optimum correlation value to perform the clustering of peaks in the identification of metabolites for NMR peaks data, which leads to the maximum
number of clusters of the variables.}
\usage{
find_corr(dataset_orig, CMETH = "pearson")
}
\arguments{
\item{dataset_orig}{List representing the dataset from an nmr peaks metabolomics experiment.}

\item{CMETH}{Correlation method used to cluster the variables. Defaults to "pearson".}
}
\value{
Value of the optimal correlation.
}
