\name{get_metabolights_study_list}
\alias{get_metabolights_study_list}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
List of metabolights studies available to download.
}
\description{
Returns a vector with the IDs of all metabolights studies available to download. Not all data in metabolights is in a format readable by specmine.
}
\usage{
get_metabolights_study_list()
}
\value{
Character vector with the IDs of the metabolights studies available to download.
}
\section{Warning}{
You must not call this function unless you have Python (>=3.5.2) installed in your machine and the module isatools.
}
\examples{
##get_metabolights_study_list()
}
