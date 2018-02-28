\name{nmr_clustering}
\alias{nmr_clustering}
\title{Cluster of variables (peaks) in a NMR peaks dataset for metabolite identification.}
\description{
Takes a dataset and performs clustering of variables, according to a correlation. The variables
will be separated into different clusters, according to a minimum correlation between variables.
Each cluster may correspond to a metabolite.
}
\usage{
nmr_clustering(dataset_orig, CMETH = "pearson", CVAL = 0.95, MVER = 2)
}
\arguments{
\item{dataset_orig}{List representing the dataset from an nmr peaks metabolomics experiment.}

\item{CMETH}{Correlation method used to cluster the variables. Defaults to "pearson".}

\item{CVAL}{Minimum correlation between variables, so they can belong to the same cluster. Can be calculated using \code{\link{find_corr}} function.}

\item{MVER}{Minimum number of variables in each cluster. Only the clusters with at least MVER variables will be returned.}
}
\value{
List with the formed clusters
}

