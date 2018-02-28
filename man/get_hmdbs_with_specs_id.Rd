\name{get_hmdbs_with_specs_id}
\alias{get_hmdbs_with_specs_id}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 HMDB code and name of HMDB 1d-NMR spectra.
}
\description{
Gives the HMDB codes and common names for the HMDB spectra ids provided.
}
\usage{
get_hmdbs_with_specs_id(spec_ids)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{spec_ids}{Vector of HMDB 1D-NMR spectra ids}
}
\value{
Dataframe with the HMDB codes and names for each HMDB spectra ids given. Each spectra ids can have more than one metabolite associated.
}
\examples{
get_hmdbs_with_specs_id(c("1022","1023", "1024", "1025"))
}
