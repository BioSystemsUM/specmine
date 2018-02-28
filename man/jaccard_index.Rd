\name{jaccard_index}
\alias{jaccard_index}
\title{Jaccard Index}
\description{
Calculates the jaccard index, i.e., the similarity between cluster and reference metabolite, in an NMR metabolite identification.
}
\usage{
jaccard_index(peaksCluster, peaksReference, PPMTOL)
}
\arguments{
    \item{peaksCluster}{
      Dataframe with the peaks of a cluster. First column corresponds to the peaks in ppm and the second one to the intensity of those peaks.
    }
    \item{peaksReference}{
      Dataframe with the peaks of a reference spectrum. First column corresponds to the peaks in ppm and the second one to the intensity of those peaks.
    }
    \item{PPMTOL}{
      ppm tolerance to consider the pairing of two peaks from the cluster and reference a match.
    }
}