\name{get_cpd_names}
\alias{get_cpd_names}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get the names of the compounds that correspond to the kegg codes given.
}
\description{
Gets the common name of the compounds of the kegg codes given.
}
\usage{
get_cpd_names(kegg_codes)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{kegg_codes}{
 Character vector with kegg codes.
}
}
\value{
Named vector with the names of the compounds. The names of the vector are the compounds' names and the vector elements the kegg codes.
}
\examples{
\donttest{
get_cpd_names(c("cpd:C00001", "cpd:C00008", "gl:G13099"))
}
}
