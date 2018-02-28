\name{create_pathway_with_reactions}
\alias{create_pathway_with_reactions}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates the pathway, with reactions included in the nodes.
}
\usage{
create_pathway_with_reactions(path, path.name, identified_cpds,
                             nodeNames="kegg", nodeTooltip=F,
                             map.zoom=F, map.layout="preset",
                             map.width=NULL, map.height=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{path}{KEGGPathway object.}
  \item{path.name}{Name of the pathway.}
  \item{identified_cpds}{Vector of kegg codes to color differently in the map.}
  \item{nodeNames}{How the nodes should be named. If "kegg", nodes are named with kegg codes. If "names", nodes are named with the common names.}
}