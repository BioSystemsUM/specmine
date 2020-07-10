\name{create_pathway_with_reactions}
\alias{create_pathway_with_reactions}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates the pathway, with reactions included in the nodes.
}
\description{
Creates a cytoscape pathway, where the reactions between compounds are also included in the nodes.
}
\usage{
create_pathway_with_reactions(path, path.name, identified_cpds,
                             nodeNames="kegg", nodeTooltip=FALSE,
                             map.zoom=FALSE, map.layout="preset",
                             map.width=NULL, map.height=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{path}{KEGGPathway object.}
  \item{path.name}{Name of the pathway.}
  \item{identified_cpds}{Vector of kegg codes to color differently in the map.}
  \item{nodeNames}{How the nodes should be named. If "kegg", nodes are named with kegg codes. If "names", nodes are named with the common names.}
  \item{nodeTooltip}{Boolean value indicating if tooltips of nodes should appear when hovering with the mouse. Does not work for all environments (e.g. can be used in shiny apps).}
  \item{map.zoom}{Boolean value indicating if a zoom widget should appear or not. Does not work for all environments (e.g. can be used in shiny apps).}
  \item{map.layout}{Layout of the map, available values are the ones of cytoscape ("breadthfirst", "preset", "cose", ...)}
  \item{map.width}{width of the map, in percentage (e.g. "80\%"). May not work as expected in some environments.}
  \item{map.height}{Height of the map, in px (e.g. "500px"). May not work as expected in some environments.}
}
