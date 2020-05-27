\name{convert_chebi_to_kegg}
\alias{convert_chebi_to_kegg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert CHEBI codes to KEGG codes.
}
\description{
Converts a vector of CHEBI codes into a vector of the corresponding KEGG codes. This is performed by using our internal library used in NMR identification, it will not have all chebi codes.
}
\usage{
convert_chebi_to_kegg(chebi_codes)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chebi_codes}{Vector with the CHEBI codes (each chebi must be structured like: CHEBI:<number>)}
}

\value{
Named vector with kegg codes and respective names. Vector names are the compound names and the vector elements the kegg codes.
}
\examples{
keggs=convert_hmdb_to_kegg(c("CHEBI:15377", "CHEBI:26078", "CHEBI:30168"))
keggs
}
