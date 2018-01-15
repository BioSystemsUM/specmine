\name{convert_hmdb_to_kegg}
\alias{convert_hmdb_to_kegg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get kegg codes from hmdb codes:
}
\usage{
convert_hmdb_to_kegg(hmdb_codes)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hmdb_codes}{
Vector with the HMDB codes (each hmdb code must have 7 digits, e.g., HMDB000001)
}
}
\value{
Named vector with kegg codes and respective names. Vector names are the compound names and the vector elements the kegg codes.
}
\examples{

}