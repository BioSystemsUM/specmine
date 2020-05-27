\name{spectra_options}
\alias{spectra_options}
\docType{data}
\title{
Information on the library of NMR reference spectra in our package.
}
\description{
This dataset provides all the information on the library of NMR spectra used as references in NMR metabolite identification.
}
\usage{data("spectra_options")}
\format{
  A data frame with 1816 observations on the following 9 variables. Each observation corresponds to a spectrum in our library.
  \describe{
    \item{\code{SPCMNS}}{a character vector with the spetcra IDs.}
    \item{\code{SPCMNM}}{a character vector with the metabolite IDs of the corresponding spectra.}
    \item{\code{FREQUENCY}}{a character vector with the frequencies under which the spectra were obtained.}
    \item{\code{NUCLEUS}}{a character vector mentioned the nucleus examined. All observations are '1H'.}
    \item{\code{PH}}{a character vector with the pH of the samples from which the spectra were obtained. May contain missing values.}
    \item{\code{TEMPERATURE}}{a character vector with the temperature under which the spectra were obtained. May contain missing values.}
    \item{\code{SOLVENT}}{a character vector with the solvent of the samples from which the spectra were obtained.}
    \item{\code{ORIGINAL_DATABASE_ID}}{whenever available, a character vector with the ID of the corresponding spectra from the database it was originally acquired from.}
    \item{\code{DATABASE}}{a character vector specifying from which database the spectra were taken from.}
  }
}
\references{
The spectra were taken from the following databases: HMDB (https://hmdb.ca), BMRB (http://www.bmrb.wisc.edu) and SDBS (https://sdbs.db.aist.go.jp). Some spectra were internally acquired and are mentioned as OUR in the DATABASE variable.
}
\examples{
data(spectra_options)
}
\keyword{datasets}
