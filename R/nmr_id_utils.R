get_kegg_groups=function(all_orgs){
  groups=unique(unlist(strsplit(all_orgs$phylogeny, ";")))
}

is_group=function(group, all_orgs){
  groups=get_kegg_groups(all_orgs)
  if (!group%in%groups) return(FALSE)
  return(TRUE)
}

is_organism=function(organism){
  if (!organism%in%names(orgs_cpds)) return(FALSE)
  return(TRUE)
}

get_organisms_in_group=function(group, all_orgs){
  if (!is_group(group, all_orgs)) stop("Invalid group.")
  
  group_search=paste(".*", group, ".*", sep="")
  organisms=as.character(all_orgs$organismCode[grep(group_search, all_orgs$phylogeny)])
  
  return(organisms)
}

compounds_in_organism=function(organism){
  
  #if (!is_organism(organism)) stop("Invalid organism code.")
  compounds=orgs_cpds[[organism]]
  
  return(compounds)
}

compounds_in_group=function(group, all_orgs){
  #Compound only needs to be in one of the organisms' group to be considered present.
  organisms=get_organisms_in_group(group, all_orgs)
  
  compounds_group=c()
  for(organism in organisms){
    compounds_organism=compounds_in_organism(organism)
    compounds_group=unique(c(compounds_group, compounds_organism))
  }
  return(compounds_group)
}

is_compound_in_organism=function(compound, organism){
  if (!is_organism(organism)) stop("Invalid organism code.")
  if (compound %in% compounds_in_organism(organism)) return(TRUE)
  return(FALSE)
}

is_compound_in_group=function(compound, group, all_orgs){
  if (!is_group(group, all_orgs)) stop("Invalid group.")
  if (compound %in% compounds_in_group(group, all_orgs)) return(TRUE)
  return(FALSE)
}

is_compound_in_entity=function(compound, entity, all_orgs){
  #entity: string of organism or group.
  #If entity is a group:
  if(is_group(entity, all_orgs)) return(is_compound_in_group(compound, entity, all_orgs))
  if(is_organism(entity)) return(is_compound_in_organism(compound, entity))
  stop("Invalid entity.")
}