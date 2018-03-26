########################################################
#####NMR METABOLITE IDENTIFICATION (FROM NMR PEAKS)#####
########################################################


#' Choose the metabolites' spectra to be the reference.
choose_nmr_references <- function(frequency, nucleus, solvent=NULL, ph=NULL, temperature=NULL){ #, org=NULL){
  ##org: organism code (or phylogeny --- later)
  
  
  env=new.env()
  data(nmr_1d_spectra, package="specmine", envir=env)
  nmr_1d_spectra=env$nmr_1d_spectra
  data(nmr_1d_spectra_options, package="specmine", envir=env)
  nmr_1d_spectra_options=env$nmr_1d_spectra_options
  
  #FREQUENCY: 400; 500; 600
  if(frequency%in%c(400, 500, 600)){
    nmr_1d_spectra=nmr_1d_spectra[nmr_1d_spectra_options$spec_id[nmr_1d_spectra_options$frequency==frequency]]
  }
  else{
    print("No frequency in library with that value. Must be 400, 500 or 600")
  }

  #NUCLEUS
  #1H or 13C
  if(nucleus%in%c("1H", "13C")){
    nmr_1d_spectra=nmr_1d_spectra[nmr_1d_spectra_options$spec_id[nmr_1d_spectra_options$nucleus==nucleus]]
  }
  else{
    print("No nucleus in library with that name. Must be '1H' or '13C'")
  }
  
  #SOLVENT
  #100%_DMSO; 5%_DMSO; acetone+DMSO+tetramethylurea; acetone+DMSO+tetramethylurea; C;
  #CCl4; CD3OD; CDCl3; cyclohexane; D2O; DMSO-d6; DMSO-d6+HCl; neat; TMS; Water
  if (!is.null(solvent)){
    if(solvent%in%levels(nmr_1d_spectra_options$solvent)){
      nmr_1d_spectra=nmr_1d_spectra[nmr_1d_spectra_options$spec_id[nmr_1d_spectra_options$solvent==solvent]]
    }
    else{
      print("No solvent in library with that name.")
    }
  }
  
  #PH
  #NA, 1.00, 1.20, 1.34, 3.00, 3.50, 4.00, 5.00, 5.50, 6.00, 6.48, 6.50, 6.90, 6.98, 6.99, 7.00
  #7.01, 7.02, 7.03, 7.10, 7.20, 7.78, 8.00, 8.50, 9.00, 10.00, 10.30, 12.00
  if(!is.null(ph)){
    if(length(ph)==2){
      nmr_1d_spectra_options=nmr_1d_spectra_options[!is.na(nmr_1d_spectra_options$ph),]
      nmr_1d_spectra=nmr_1d_spectra[nmr_1d_spectra_options$spec_id[nmr_1d_spectra_options$ph<=ph[2]&nmr_1d_spectra_options$ph>=ph[1]]]
    }
    else if (length(ph)==1){
      nmr_1d_spectra=nmr_1d_spectra[nmr_1d_spectra_options$spec_id[nmr_1d_spectra_options$ph==ph]]
    }
  }
  
  #TEMPERATURE
  #25, 50, NA
  if (!is.null(temperature)){
    if(temperature%in%c(25,50)){
      nmr_1d_spectra=nmr_1d_spectra[nmr_1d_spectra_options$spec_id[nmr_1d_spectra_options$temperature==temperature]]
    }
    else{
      print("No temperature in library with that value.")
    }
  }
  return(nmr_1d_spectra[!is.na(names(nmr_1d_spectra))])
}



#' GIVE THE METABOLITES HMDBS and names TO WHICH THE SPECTRA IDS BELONG
get_hmdbs_with_specs_id=function(spec_ids){
  env=new.env()
  data(conversion_table, package="specmine", envir=env)
  hmdbs_with_spec_refs=env$conversion_table[!is.na(env$conversion_table$SPECTRA_NMR_ONED_OWN), c("HMDB", "NAME", "SPECTRA_NMR_ONED_OWN")]
  hmdbs=c()
  spec=c()
  names=c()
  for(spec_id in spec_ids){
    for(i_spec in 1:length(hmdbs_with_spec_refs$SPECTRA_NMR_ONED_OWN)){
      if(length(grep(paste(".*", spec_id, ".*", sep=""), hmdbs_with_spec_refs$SPECTRA_NMR_ONED_OWN[i_spec]))>0){
        hmdbs=c(hmdbs, as.character(hmdbs_with_spec_refs$HMDB[i_spec]))
        names=c(names, as.character(hmdbs_with_spec_refs$NAME[i_spec]))
        spec=c(spec, spec_id)
      }
    }
  }
  res=data.frame(hmdbs, names, spec)
  return(res)
}




#' FIND THE BEST CORRELATION VALUE FOR THE CLUSTERING OF PEAKS
find_corr <- function(dataset_orig, CMETH='pearson', maxPeaks=40) {
  
  #CODE ADAPTED FROM THE CODE USED IN THE ARTICLE "An efficient spectra processing method for
  #metabolite identification from 1H-NMR metabolomics data", by Daniel Jacob, Catherine Deborde,
  #Annick Moing. 
  
  #Dataset normalization:
  dataset_norm=normalize(dataset_orig, "median")
  matrix=t(dataset_norm$data)
  colnames(matrix)=paste("V",colnames(matrix), sep="")
  
  #Correlation matrix between variables:
  cor_mat <- cor(matrix,method=CMETH)
  cor_mat[ lower.tri(cor_mat, diag=TRUE) ]<- 0
  
  CVAL_MIN<-0.9
  CVAL_MAX<-0.9999
  CVAL_STEP<-0.0001
  CVAL_OPT_MIN<-0
  CVAL_OPT_MAX<-0
  CVAL_OPT<-0
  NBCLUSTERS_MAX<-0
  CVAL_CRIT<-0
  
  MAXSIZE<-maxPeaks
  
  cat("CVAL","nb_clusters","nb_clusters_2","size_max","Criterion","nb_buckets",sep=";"); cat ("\n")
  for (CVAL in seq(CVAL_MIN, CVAL_MAX, by=CVAL_STEP)) {
    cor_mati <- cor_mat
    cor_mati[ cor_mati < CVAL] <- 0
    graph <- igraph::graph.adjacency(cor_mati>CVAL, weighted=TRUE, mode="upper")
    igraph::E(graph)$weight<-t(cor_mati)[t(cor_mati)>CVAL]
    igraph::V(graph)$label<- igraph::V(graph)$name
    cliques <- sapply(igraph::decompose.graph(graph), igraph::vcount)
    ordcli <- order(cliques,decreasing = T)
    nb_clusters<-0
    nb_clusters_2<-0
    nb_buckets<-0
    for (i in 1:length(ordcli)) {
      if (cliques[ordcli[i]]>=2) {
        nb_clusters <- nb_clusters + 1
        nb_buckets <- nb_buckets + cliques[ordcli[i]]
      }
      if (cliques[ordcli[i]]==2)
        nb_clusters_2 <- nb_clusters_2 + 1
    }
    size_max <- cliques[ordcli[1]]
    CRIT<-size_max/nb_clusters
    cat(CVAL,nb_clusters,nb_clusters_2,size_max,-20*log10(CRIT),nb_buckets,sep=";"); cat ("\n")
    if (size_max>MAXSIZE) next
    if (size_max<=MAXSIZE && CVAL_OPT_MIN==0) {
      CVAL_OPT_MIN<-CVAL
      CVAL_CRIT<-CRIT
      CVAL_OPT_MAX<-CVAL
      CVAL_OPT<-CVAL
      NBCLUSTERS_MAX<-size_max
      next
    }
    if (CRIT<CVAL_CRIT) {
      CVAL_CRIT<-CRIT
      CVAL_OPT_MAX<-CVAL
      if (size_max>NBCLUSTERS_MAX) {
        NBCLUSTERS_MAX<-size_max
        CVAL_OPT<-CVAL
      }
      next
    }  
  }
  cat("# -----------------\n")
  cat("# Optimum Corr Min=",CVAL_OPT_MIN,", Max=",CVAL_OPT_MAX,", Mean=",CVAL_OPT,"\n",sep="")
  
  return(CVAL_OPT)
}





#' CLUSTER OF VARIABLES (PEAKS)
nmr_clustering <- function(dataset_orig, CMETH='pearson', CVAL=0.95, MVER=2)
{
  #CODE ADAPTED FROM THE CODE USED IN THE ARTICLE "An efficient spectra processing method for
  #metabolite identification from 1H-NMR metabolomics data", by Daniel Jacob, Catherine Deborde,
  #Annick Moing. 
  
  final.results=list()
  
  dataset_norm=normalize(dataset_orig, "median")
  matrix=t(dataset_norm$data)
  colnames(matrix)=paste("V",colnames(matrix), sep="")
  
  #Correlation matrix between variables:
  cor_mat <- cor(matrix,method=CMETH)
  cor_mat[lower.tri(cor_mat, diag=TRUE)]<- 0
  cor_mat[cor_mat < CVAL] <- 0
  
  #vertices: peaks, edges: links two edges with a correlation value bigger than CVAL
  graph <- igraph::graph.adjacency(cor_mat>CVAL, weighted=TRUE, mode="upper")
  #weight of edges will be the correlation between them
  igraph::E(graph)$weight<-t(cor_mat)[t(cor_mat)>CVAL]
  #Label will be the names of the vertices
  igraph::V(graph)$label<- igraph::V(graph)$name
  #Separates graph in the different sub graphs (a sub graph is not connected with anny other
  #subgraph - no correlation) and counts the number of peaks in each subgraph
  cliques <- sapply(igraph::decompose.graph(graph), igraph::vcount)
  ordcli <- order(cliques,decreasing = T) #Dá os indices dos grafos, ordenados pelo grafo com
                                        #maior nº de vertices para o menor
  M <- NULL
  g<-0
  
  data_orig=t(dataset_orig$data)
  colnames(data_orig)=paste("V",colnames(data_orig), sep="")
  for (i in 1:length(ordcli)) {
    ind<-ordcli[i]
    if (cliques[ind]>=MVER) {
      g <- g + 1
      subg <- igraph::decompose.graph(graph)[[ind]]
      M <- data_orig[,colnames(data_orig) %in% igraph::V(subg)$name]
      clust <- cbind(as.vector(as.numeric(lapply(igraph::V(subg)$name, function(x) substr(x,2,nchar(x))))),
                     as.vector(apply(t(M),1,mean)))
      final.results[[paste("Cluster", g)]]=clust
    }
  }
  
  return(final.results)
}





#' JACCARD INDEX
#' 
#' Calculates the jaccard index, i.e., the similarity between cluster and reference metabolite.
#' 
jaccard_index=function(peaksCluster, peaksReference, PPMTOL){
  
  res=list()
  
  #Save the peaks from the cluster and references that matched:
  matched_peaks_clust=c()
  matched_peaks_ref=c()
  
  matched_peaks=0
  
  #Cluster and reference peaks:
  clust_cs=sort(peaksCluster[,1])
  ref_cs=sort(peaksReference$chemical.shift)
  
  #For each cluster peak, compare it to the reference peaks:
  last_pR=1
  for(pC in 1:length(clust_cs)){
    if (last_pR>length(ref_cs)) break
    for(pR in last_pR:length(ref_cs)){
      inf_limit=ref_cs[pR]-PPMTOL
      if( round(clust_cs[pC],3) < round(inf_limit,3)) break
      sup_limit=ref_cs[pR]+PPMTOL
      if( round(clust_cs[pC],3) <= round(sup_limit,3)){
        matched_peaks = matched_peaks + 1
        last_pR=pR+1
        
        matched_peaks_clust=c(matched_peaks_clust, round(clust_cs[pC],3))
        matched_peaks_ref=c(matched_peaks_ref, round(ref_cs[pR],3))
        break
      }
    }
  }
  score = matched_peaks / (length(unique(ref_cs))+length(unique(clust_cs))-matched_peaks)
  
  #Store results:
  res$score=score
  res$matched_peaks_ref=matched_peaks_ref
  res$matched_peaks_clust=matched_peaks_clust
  res$reference_peaks=ref_cs
  
  
  return(res)
}





#' FUNCTION TO DO THE NMR IDENTIFICATION
#' 
#' This function performs metabolite identification on a dataset of nmr peaks, finding the top
#' 5 reference metabolites that matched with each cluster of variables formed.
nmr_identification <- function(dataset, ppm.tol=0.03,
                               clust.method='pearson', clust.treshold=NULL, clust.peaks.min=2,
                               clust.maxPeaks=40, clust.nTop=5,
                               freq=500, nucl="1H", solv=NULL, pH=NULL, temp=NULL){
  
  #Choose reference metabolites:
  cat("Getting reference metabolites...\n")
  references=choose_nmr_references(frequency=freq, nucleus=nucl, solvent=solv, ph=pH,
                                   temperature=temp)
  #Get names of the reference metabolites:
  hmdbs_to_spec=get_hmdbs_with_specs_id(names(references))
  
  #Find best correlation value for the clusters correlation treshold, if a value is not given:
  if (is.null(clust.treshold)){
    cat("Getting best correlation value...\n")
    if (is.null(clust.maxPeaks)){
      nPeaks=c()
      for(spectra in references){
        nPeaks=c(nPeaks, dim(spectra)[1])
      }
      nPeaks=max(nPeaks)
    }
    else nPeaks=clust.maxPeaks
    clust.treshold=find_corr(dataset, CMETH=clust.method, maxPeaks=nPeaks)
  }
  #Get clusters:
  cat("Getting the clusters of the peaks...\n")
  clusters=nmr_clustering(dataset, CMETH=clust.method, CVAL=clust.treshold,
                                  MVER=clust.peaks.min)
  
  #Results variable:
  res=list()
  
  #For each cluster, compare with each reference:
  cat("Matching clusters with reference metabolites...\n")
  for (clust_id in 1:length(clusters)){
    full_res=list()
    score_clust=c()
    n_score_clust=c()
    for (ref_id in 1:length(references)){
      res_ref=jaccard_index(clusters[[clust_id]], references[[ref_id]], PPMTOL=ppm.tol)
      
      spec_id=names(references)[ref_id]
      hmdb_id=as.character(hmdbs_to_spec$hmdbs[hmdbs_to_spec$spec==spec_id])
      for(h in hmdb_id){
        full_res[[h]]=res_ref
        score_clust=c(score_clust, res_ref$score)
        n_score_clust=c(n_score_clust, h)
      }
    }
    #Store only the top clust.nTop matches of the cluster and the cluster peaks
    names(score_clust)=n_score_clust
    score_clust=sort(score_clust, decreasing=T)[1:clust.nTop]
    for (i in 1:length(score_clust)){
      if (score_clust[i]==0){
        take=i:length(score_clust)
        score_clust=score_clust[-take]
        break
      }
    }
    full_res=full_res[names(score_clust)]
    
    res[[paste("Cluster", clust_id, sep="")]]$cluster.peaks=clusters[[clust_id]]
    res[[paste("Cluster", clust_id, sep="")]]$metabolites.matched=full_res
    res[[paste("Cluster", clust_id, sep="")]]$summary=score_clust
  }
  cat("Done.\n")
  
  return(res)
}


