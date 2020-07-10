#####################
#####CONVERSIONS#####
#####################

#' Get the names of the compounds that correspond to the kegg codes given:
get_cpd_names=function(kegg_codes){
  if(!requireNamespace("KEGGREST", quietly = TRUE)){
    stop("Package KEGGREST needed for this function to work. Please install it: BiocManager::install('KEGGREST').",
         call. = FALSE)
  }
  
  n=c()
  n_kegg=c()
  cpds_to_get=c()
  #Get names from the ones that are present in the conversion table:
  for (i in 1:length(kegg_codes)){
    kegg=kegg_codes[i]
    tab_kegg=strsplit(kegg, ":")[[1]][2]
    for(i in 1:length(codes$KEGG)){
      kegg_ref=codes$KEGG[i]
      if (tab_kegg %in% unlist(strsplit(kegg_ref, "; ")[[1]])){
        n=c(n,codes$NAME[i])
        n_kegg=c(n_kegg, kegg)
      }
      else{
        cpds_to_get=c(cpds_to_get, kegg)
      }
    }
  }
  
  #If not present in conversion table, get from kegg itself:
  if(!is.null(cpds_to_get)){
    all_cpds=c(KEGGREST::keggList("compound"), KEGGREST::keggList("glycan"))
    cpds_names_kegg=all_cpds[cpds_to_get]
    all_n_kegg=strsplit(cpds_names_kegg, "; ")
    for (i in 1:length(all_n_kegg)){
      kegg_name=all_n_kegg[[i]][1]
      kegg_code=names(all_n_kegg)[i]
      n=c(n, kegg_name)
      n_kegg=c(n_kegg, kegg_code)
    }
  }
  
  names(n_kegg)=n
  return(n_kegg)
}

#' Get kegg codes from hmdb codes:
convert_hmdb_to_kegg=function(hmdb_codes){
  
  hmdbs=toupper(hmdb_codes)
  names_cpds=c()
  keggs=c()
  for (hmdb in hmdbs){
    for(i in 1:length(codes$HMDB)){
      hmdb_ref=codes$HMDB[i]
      if (hmdb %in% unlist(strsplit(hmdb_ref, "; ")[[1]])){
        keggs_for_hmdb = strsplit(codes$KEGG[i], "; ")[[1]]
        if(length(keggs_for_hmdb)!=0){
          cpd=paste("cpd:", keggs_for_hmdb, sep="")
          keggs=c(keggs, cpd)
          names_cpds=c(names_cpds, rep(codes$NAME[codes$HMDB==hmdb], length(cpd)))
        }
      }
    }
  }
  names(keggs)=names_cpds
  return(keggs)
}

#' Get kegg codes from chebi codes:
convert_chebi_to_kegg=function(chebi_codes){
  
  chebis=toupper(chebi_codes)
  names_cpds=c()
  keggs=c()
  for (chebi in chebis){
    for(i in 1:length(codes$HMDB)){
      chebi_ref=codes$CHEBI[i]
      if (chebi %in% unlist(strsplit(chebi_ref, "; ")[[1]])){
        keggs_for_chebi = strsplit(codes$KEGG[i], "; ")[[1]]
        if(length(keggs_for_chebi)!=0){
          cpd=paste("cpd:", keggs_for_chebi, sep="")
          keggs=c(keggs, cpd)
          names_cpds=c(names_cpds, rep(codes$NAME[codes$CHEBI==chebi], length(cpd)))
        }
      }
    }
  }
  names(keggs)=names_cpds
  return(keggs)
}

#' Get kegg codes from spcmnm codes:
convert_multiple_spcmnm_to_kegg=function(spcmnm_codes){
  
  spcmnms=toupper(spcmnm_codes)
  names_cpds=c()
  keggs=c()
  for (spcmnm in spcmnms){
    if(!spcmnm%in%codes$SPCMNM) stop("The following SPCMNM code is not available: ", spcmnm, "\n")
    keggs_for_spcmnm=unlist(strsplit(codes$KEGG[codes$SPCMNM==spcmnm], "; ")[[1]])
    if(length(keggs_for_spcmnm)!=0){
      cpd=paste("cpd:", keggs_for_spcmnm, sep="")
      keggs=c(keggs, cpd)
      names_cpds=c(names_cpds, rep(codes$NAME[codes$SPCMNM==spcmnm], length(cpd)))
    }
  }
  names(keggs)=names_cpds
  return(keggs)
}


###########################
#####ORGANISMS - PATHS#####
###########################

#' Get code, t number, full name and phylogeny of all organisms in KEGG:
get_OrganismsCodes=function(){
  if(!requireNamespace("KEGGREST", quietly = TRUE)){
    stop("Package KEGGREST needed for this function to work. Please install it: BiocManager::install('KEGGREST').",
         call. = FALSE)
  }
  organisms=KEGGREST::keggList("organism")
  organisms.Inf=data.frame(Tnumber=organisms[,1], organismCode=organisms[,2],
                           speciesNames=organisms[,3], phylogeny=organisms[,4], stringsAsFactors=FALSE)
  return(organisms.Inf)
}

#' Get vector with paths numbers that occur in the given organism, named with the full path name:
get_metabPaths_org=function(org_code){
  if(!requireNamespace("KEGGREST", quietly = TRUE)){
    stop("Package KEGGREST needed for this function to work. Please install it: BiocManager::install('KEGGREST').",
         call. = FALSE)
  }
  
  if (!org_code%in%get_OrganismsCodes()[["organismCode"]]) stop("Invalid organism code.")
  
  paths=unique(KEGGREST::keggLink("pathway", org_code))
  paths_n=substr(paths, nchar(paths)-5+1, nchar(paths))
  
  maps_take_out=c("00190", "00195", "00196", "01100", "01110", "01120", "01130", "01200", "01210", "01212", "01230", "01220", "00511", "00514",
                 "00533", "01010", "01060", "01062", "01063", "01064", "01065", "01066", "01070", "01052", "01054", "00121")
  idx_take_out=as.vector(na.exclude(match(maps_take_out, paths_n)))
  idx_take_out=c(idx_take_out, as.vector(na.exclude(match(paths_n[as.numeric(paths_n)>1230], paths_n))))
  paths_n=paths_n[-idx_take_out]
  
  pathways_all=KEGGREST::keggList("pathway")
  paths_names=pathways_all[match(paste("path:map", paths_n, sep=""), names(pathways_all))]
  paths_n=paste(org_code, paths_n, sep="")
  names(paths_n)=paths_names
  
  return(paths_n)
  
}



#######################################
#####ORGANISMS - PATHS - COMPOUNDS#####
#######################################

#' Get only the paths of the organism that contain given compounds:
get_paths_with_cpds_org=function(organism_code, compounds, full.result=TRUE){
  if(!requireNamespace("KEGGgraph", quietly = TRUE)){
    stop("Package KEGGgraph needed for this function to work. Please install it: BiocManager::install('KEGGgraph').",
         call. = FALSE)
  }
  
  org_paths=get_metabPaths_org(organism_code)
  
  paths_with_cpds=c()
  names_paths=c()
  compounds_in_paths=c()
  compounds_in_paths_names=c()
  ratio=c()
  for (i in 1:length(org_paths)){
    path=org_paths[i]
    path_name=names(org_paths)[i]
    #message("Analysing pathway", path, "\n")
    path_r=get_MetabolitePath(path)
    nods=KEGGgraph::nodes(convert_keggpathway_2_reactiongraph(path_r))
    cpds_idx=as.vector(na.exclude(match(nods, compounds)))
    if (length(cpds_idx)>0){
      paths_with_cpds=c(paths_with_cpds, path)
      names_paths=c(names_paths, path_name)
      ratio=c(ratio, length(cpds_idx)/length(nods))
      if (full.result){
        compounds_in_paths=c(compounds_in_paths, paste(compounds[unique(cpds_idx)], collapse="; "))
        compounds_in_paths_names=c(compounds_in_paths_names, paste(names(compounds)[unique(cpds_idx)], collapse="; "))
      }
    }
  }
  order_idx=order(ratio)
  if(full.result) return(data.frame(pathways=paths_with_cpds, ratio=ratio, compounds=compounds_in_paths, compounds_names=compounds_in_paths_names, row.names=names_paths, stringsAsFactors=FALSE)[order_idx,])
  return(data.frame(pathways=paths_with_cpds, ratio=ratio, row.names=names_paths, stringsAsFactors=FALSE)[order_idx,])
  
}



########################
#####CREATE PATHWAY#####
########################

#' Returns an object of KEGGPathway of the pathway especified in pathcode
get_MetabolitePath=function(pathcode){
  pathObj=KEGGREST::keggGet(pathcode, option="kgml")
  parsedPath=KEGGgraph::parseKGML(pathObj)
  return(parsedPath)
}

#' Convert KEGGPathway object to graph object
convert_keggpathway_2_reactiongraph=function(pathObj){
  reactionGraphObj=KEGGgraph::KEGGpathway2reactionGraph(pathObj)
  return(reactionGraphObj)
}

#' Creates the pathway, with reactions included in the nodes
create_pathway_with_reactions=function(path, path.name, identified_cpds,
                                       nodeNames="kegg", nodeTooltip=FALSE,
                                       map.zoom=FALSE, map.layout="preset",
                                       map.width=NULL, map.height=NULL){
  
  #Vectors to construct edges:
  source=c()
  target=c()
  edgeSourceShape=c()
  edgeTargetShape=c()
  colorLine=c()
  
  #Vectors to construct nodes:
  nodes=c()
  shape=c()
  height=c()
  width=c()
  color=c()
  
  #Get nodes and edges of reactions:
  reaction_map=KEGGgraph::getReactions(path)
  
  for (reaction in reaction_map){ #For each reaction in path(s)
    
    n=KEGGgraph::getName(reaction) #Get reaction name
    names=strsplit(n, " ")[[1]]
    for (name in names){
      
      t=KEGGgraph::getType(reaction) #Reaction type: reversible or irreversible
      substrate=c()
      for (l in strsplit(KEGGgraph::getSubstrate(reaction), " ")) substrate=c(substrate, l) #Get substrates
      n_substrate=length(substrate)
      product=c()
      for (l in strsplit(KEGGgraph::getProduct(reaction), " ")) product=c(product, l) #Get products
      n_product=length(product)
      n_compounds=n_substrate+n_product
      n_edges=n_compounds+1
      
      source=c(source, substrate)
      target=c(target, rep(name, n_substrate))
      source=c(source, rep(name, n_product))
      target=c(target, product)
      
      #If reaction is irreversible, edge line will be red:
      if(t=="irreversible"){
        colorLine=c(colorLine, rep("#bf000c", n_compounds))
        color=c(color, rep("#888888", n_compounds), "#bf000c")
      } 
      else{ #If reaction is reversible, edge line will be green:
        colorLine=c(colorLine, rep("#145b02", n_compounds))
        color=c(color, rep("#888888", n_compounds), "#145b02")
      }
      #edgeSourceShape=c(edgeSourceShape, rep("none", n_compounds))
      #edgeTargetShape=c(edgeTargetShape, rep("triangle", n_compounds))
      
      nodes=c(nodes, substrate, product, name)
      #ellipse for compound nodes and roundrectangle for reaction nodes:
      shape=c(shape, rep("ellipse", n_compounds), "roundrectangle")
      height=c(height, rep(80, n_compounds), 30)
      width=c(width, rep(80, n_compounds), 80)
    }
  }
  
  #Get possible maps to which the present map conects to, and add it to nodes and reactions:
  map_names=c()
  maps=c()
  path_p=substr(path.name, nchar(path.name)-5+1, nchar(path.name))
  maps_con=maps_con[grep(paste(".*", path_p, ".*", sep=""), maps_con$in_map),]
  for (node in path@nodes){
    if(KEGGgraph::getType(node)=="map" & !(KEGGgraph::getName(node)[1]%in%maps) & length(grep("^TITLE:", node@graphics@name))==0){
      map_name=KEGGgraph::getName(node)
      map_names=c(map_names, node@graphics@name)
      maps=c(maps, map_name)
      
      path_n=substr(map_name, nchar(map_name)-4, nchar(map_name))
      
      cpds=maps_con$node[maps_con$map==path_n]
      cpds_ant=c()
      if (sum(cpds%in%nodes)!=0){
        #graph information for path:
        nodes=c(nodes, map_name)
        color=c(color, "#D2691E")
        shape=c(shape, "rectangle")
        height=c(height, 80)
        width=c(width, 80)
        for(cpd in cpds){
          if(cpd%in%nodes & !(cpd%in%cpds_ant)){
            type=maps_con$type[maps_con$map==path_n & maps_con$node==cpd]
            cpds_ant=c(cpds_ant, cpd)
            if(type%in%c("reversible", "in")){
              #graph information for other node:
              nodes=c(nodes, cpd)
              color=c(color, "#888888")
              shape=c(shape, "ellipse")
              height=c(height, 80)
              width=c(width, 80)
              #graph information for edge:
              source=c(source, cpd)
              target=c(target, map_name)
              if(type=="reversible") colorLine=c(colorLine, "#145b02")
              else colorLine=c(colorLine, "#bf000c")
            }
            else{
              #graph information for other node:
              nodes=c(nodes, cpd)
              color=c(color, "#888888")
              shape=c(shape, "ellipse")
              height=c(height, 80)
              width=c(width, 80)
              #graph information for edge:
              source=c(source, map_name)
              target=c(target, cpd)
              colorLine=c(colorLine, "#bf000c")
            }
          }
        }
      }
    }
  }
  names(map_names)=maps
  
  #Names of the nodes
  if(nodeNames=="kegg"){
    name=c()
    new_name=strsplit(nodes, ":")
    for (n in new_name){name=c(name, n[2])}
  }
  else if(nodeNames=="names"){
    name=nodes
    cpds=grep("^cpd:", name, value=TRUE)
    cpds=c(cpds, grep("^gl:", name, value=TRUE))
    not_cpds=grep("^cpd:", name, value=TRUE, invert=TRUE)
    named_cpds=get_cpd_names(cpds)
    for (cpd in named_cpds) name[name==cpd]=names(named_cpds)[named_cpds==cpd]
    for (code in not_cpds){
      if(length(grep("^path:", code, value=TRUE))) name[name==code]=map_names[code]
      name[name==code]=strsplit(code, ":")[[1]][2]
    }
  }
  
  #Set positions from path if preset layout is chosen:
  if (map.layout=="preset"){
    node_name=c()
    x=c()
    y=c()
    x.data=c()
    y.data=c()
    
    for (node in path@nodes){
      x=c(x, node@graphics@x)
      y=c(y, node@graphics@y)
      if(is.na(node@reaction)) node_name=c(node_name, paste(node@name, collapse=" "))
      else if(KEGGgraph::getType(node)=="map") node_name=c(node_name, paste(node@name, collapse=" "))
      else node_name=c(node_name, node@reaction)
    }
    position.df=data.frame(x=x, y=y, node_name=node_name)
    position.df$x=position.df$x*2
    position.df$y=position.df$y*2
    
    for (node in nodes){
      x=position.df$x[position.df$node_name==node][1]
      y=position.df$y[position.df$node_name==node][1]
      
      if(sum(is.na(c(x,y)))>0){
        target_x=position.df$x[position.df$node_name==target[match(node, source)]][1]
        target_y=position.df$y[position.df$node_name==target[match(node, source)]][1]
        source_x=position.df$x[position.df$node_name==source[match(node, target)]][1]
        source_y=position.df$y[position.df$node_name==source[match(node, target)]][1]
        x=(target_x+source_x)/2
        y=(target_y+source_y)/2
      }
      x.data=c(x.data, x)
      y.data=c(y.data, y)
    }
    
    #Create node and edge data:
    nodeData=data.frame(id=nodes, name=name, shape, height=height, width=width, color, nodeLabelColor=rep("#000000", length(nodes)), x=x.data, y=y.data, stringsAsFactors=FALSE)
    edgeData=data.frame(source, target, color=colorLine, label=rep("", length(source)), labelColor=rep("#888888", length(source)), stringsAsFactors=FALSE)#, edgeSourceShape, edgeTargetShape)
  }
  else{
    #Create node and edge data:
    nodeData=data.frame(id=nodes, name=name, shape, height=height, width=width, color, nodeLabelColor=rep("#000000", length(nodes)), stringsAsFactors=FALSE)
    edgeData=data.frame(source, target, color=colorLine, label=rep("", length(source)), labelColor=rep("#888888", length(source)), stringsAsFactors=FALSE)#, edgeSourceShape, edgeTargetShape)
  }
  
  #Change color nodes of the given compounds:
  cpds_to_color=c()
  for (cpd in identified_cpds){
    if (cpd%in%nodes) cpds_to_color=c(cpds_to_color, cpd)
  }
  for (cpd in cpds_to_color){ #color the identified compound nodes
    nodeData$color[nodeData$id==cpd]="#0061ff"
  }
  
  
  #Create the map:
  network <- rcytoscapejs::createCytoscapeJsNetwork(nodeData, edgeData)
  rcytoscapejs::rcytoscapejs(network$nodes, network$edges, showPanzoom=map.zoom, layout=map.layout, width=map.width, height=map.height)#, height="2000px")
}




#' Creates the pathway wanted. If any of the given compounds is present in the pathway, it is coloured differently. 
pathway_analysis=function(compounds, pathway, #cpd.type="kegg",
                          nodeNames="kegg", nodeTooltip=FALSE,
                          map.zoom=FALSE, map.layout="preset", map.width=NULL, map.height=NULL){
  #Path code: "hsa00010", for example
  #For now, "compounds" has to be in kegg codes (cpd.type="kegg").
  if (!requireNamespace("KEGGgraph", quietly = TRUE)) {
    if(!requireNamespace("KEGGREST", quietly = TRUE)){
      if(!requireNamespace("rcytoscapejs", quietly = TRUE)){
        stop("Packages KEGGgraph, KEGGREST and rcytoscapejs are needed for this function to work. Please install them: BiocManager::install(c('KEGGgraph','KEGGREST')). For rcytoscape, use: devtools::install_github('cytoscape/r-cytoscape.js@v0.0.7')",
             call. = FALSE)
      }
      stop("Packages KEGGgraph and KEGGREST are needed for this function to work. Please install them: BiocManager::install(c('KEGGgraph','KEGGREST')).",
           call. = FALSE)
    }
    else if(!requireNamespace("rcytoscapejs", quietly = TRUE)){
      stop("Packages KEGGgraph and rcytoscapejs are needed for this function to work. Please install them: BiocManager::install(c('KEGGgraph')). For rcytoscape, use: devtools::install_github('cytoscape/r-cytoscape.js@v0.0.7')",
           call. = FALSE)
    }
    stop("Package KEGGgraph needed for this function to work. Please install it: BiocManager::install('KEGGgraph').", call. = FALSE)
  }
  else if(!requireNamespace("KEGGREST", quietly = TRUE)){
    if(!requireNamespace("rcytoscapejs", quietly = TRUE)){
      stop("Packages KEGGREST and rcytoscapejs are needed for this function to work. Please install them: BiocManager::install(c('KEGGREST')). For rcytoscape, use: devtools::install_github('cytoscape/r-cytoscape.js@v0.0.7')",
           call. = FALSE)
    }
    stop("Package KEGGREST needed for this function to work. Please install it: BiocManager::install('KEGGREST').",
         call. = FALSE)
  }
  
  #Get pathway maps and graphs
  message("Getting pathway map\n")
  pathMap=get_MetabolitePath(pathway)
  reactionObj=convert_keggpathway_2_reactiongraph(pathMap)
  
  #Create the pathway map:
  message("Creating pathway\n")
  create_pathway_with_reactions(pathMap, pathway, compounds,
                                nodeNames, nodeTooltip,
                                map.zoom, map.layout, map.width, map.height)
  
}
