\name{convert_multiple_spcmnm_to_kegg}
\alias{convert_multiple_spcmnm_to_kegg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert specmine metabolite codes to KEGG codes.
}
\description{
Converts a vector of specmine metabolite codes into a vector of the corresponding KEGG codes. This is performed by using our internal library used in NMR identification, it will not have all chebi codes.
}
\usage{
convert_multiple_spcmnm_to_kegg(spcmnm_codes)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{spcmnm_codes}{Vector with the SPCMNM codes (each chebi must be structured like: SPCMNM<number>)}
}
\value{
Named vector with kegg codes and respective names. Vector names are the compound names and the vector elements the kegg codes.
}
\examples{
keggs=convert_multiple_spcmnm_to_kegg(c("SPCMNM2111", "SPCMNM2142", "SPCMNM069774"))
keggs
}