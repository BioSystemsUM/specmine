\name{convert_hmdb_to_kegg}
\alias{convert_hmdb_to_kegg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert HMDB codes to KEGG codes.
}
\description{
Converts a vector of HMDB codes into a vector of the corresponding KEGG codes. This is performed by using our internal library used in NMR identification, it will not have all chebi codes.
}
\usage{
convert_hmdb_to_kegg(hmdb_codes)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hmdb_codes}{
Vector with the HMDB codes (each hmdb code must have 7 digits, e.g., HMDB0000001)
}
}
\value{
Named vector with kegg codes and respective names. Vector names are the compound names and the vector elements the kegg codes.
}
\examples{
keggs=convert_hmdb_to_kegg(c("HMDB0000001", "HMDB0000008", "HMDB0000246"))
keggs
}
