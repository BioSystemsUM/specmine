########################################################
#####NMR METABOLITE IDENTIFICATION (FROM NMR PEAKS)#####
########################################################


#' CHOOSE THE METABOLITES TO BE THE REFERENCE
#' 
#' This function returns the reference spectra, according to the characteristics chosen
#' (frequency, nucleus, solvent, pH and temperature). Only the frequency and nucleus arguments
#' are obligatory.
#' 
#' @param frequency Frequency of the reference spectra, in Hz. Must either be 15, 22, 25, 50, 90, 100, 300, 400, 500 or 600.
#' @param nucleus Atomic nuclei. Possible values: "1H", "13C"
#' @param solvent If given, the solvent of the sample. Possbile values: "100\%_DMSO"; "5\%_DMSO"; "acetone+DMSO+tetramethylurea"; "acetone+DMSO+tetramethylurea"; "C"; "CCl4"; "CD3OD"; "CDCl3"; "cyclohexane"; "D2O"; "DMSO-d6"; "DMSO-d6+HCl"; "neat"; "TMS"; "Water".
#' @param ph If given, a number corresponding to the sample's pH or a vector of length two indicating a pH interval.
#' @param temperature If given, the temperature of the sample, in Celsius. Must either be 25 or 50.
#' 
#' @return Returns a list with the reference spectra that have the chosen characteristics. Each spectra is a list with two elements: a list of the chemical shifts and and a list of the respective intensities.
#' 
#' @examples 
#' ##Gives all the reference spectra with a frequency of 500 Hz and the nucleus 1H:
#' ref=choose_nmr_references(500, "1H")
#' 
#' @export
choose_nmr_references <- function(frequency, nucleus, solvent=NULL, ph=NULL, temperature=NULL){
  
  data(nmr_1d_spectra, package="specmine")
  data(nmr_1d_spectra_options, package="specmine")
  
  #FREQUENCY
  #15 (15.08; 15.09); 22 (22.5; 22.53); 25 (25.16); 50 (50.18; 50.32); 90;
  #100 (100.40; 100.41; 100.54; 100.7; 100.72); 300; 400; 500; 600
  if(frequency==15){
    nmr_1d_spectra=nmr_1d_spectra[rownames(nmr_1d_spectra_options)[nmr_1d_spectra_options$frequency<16]]
  }
  else if(frequency==22){
    nmr_1d_spectra=nmr_1d_spectra[rownames(nmr_1d_spectra_options)[nmr_1d_spectra_options$frequency<23&nmr_1d_spectra_options$frequency>=22]]
  }
  else if(frequency==25){
    nmr_1d_spectra=nmr_1d_spectra[rownames(nmr_1d_spectra_options)[nmr_1d_spectra_options$frequency<26&nmr_1d_spectra_options$frequency>=25]]
  }
  else if(frequency==50){
    nmr_1d_spectra=nmr_1d_spectra[rownames(nmr_1d_spectra_options)[nmr_1d_spectra_options$frequency<51&nmr_1d_spectra_options$frequency>=50]]
  }
  else if(frequency==100){
    nmr_1d_spectra=nmr_1d_spectra[rownames(nmr_1d_spectra_options)[nmr_1d_spectra_options$frequency<101&nmr_1d_spectra_options$frequency>=100]]
  }
  else if(frequency%in%c(90, 300, 400, 500, 600)){
    nmr_1d_spectra=nmr_1d_spectra[rownames(nmr_1d_spectra_options)[nmr_1d_spectra_options$frequency==frequency]]
  }
  else{
    print("No frequency in library with that value. Must be 15, 22, 25, 50, 90, 100, 300, 400, 500 or 600")
  }

  #NUCLEUS
  #1H or 13C
  if(nucleus%in%c("1H", "13C")){
    nmr_1d_spectra=nmr_1d_spectra[rownames(nmr_1d_spectra_options)[nmr_1d_spectra_options$nucleus==nucleus]]
  }
  else{
    print("No nucleus in library with that name. Must be '1H' or '13C'")
  }
  
  #SOLVENT
  #100%_DMSO; 5%_DMSO; acetone+DMSO+tetramethylurea; acetone+DMSO+tetramethylurea; C;
  #CCl4; CD3OD; CDCl3; cyclohexane; D2O; DMSO-d6; DMSO-d6+HCl; neat; TMS; Water
  if (!is.null(solvent)){
    if(solvent%in%levels(nmr_1d_spectra_options$solvent)){
      nmr_1d_spectra=nmr_1d_spectra[rownames(nmr_1d_spectra_options)[nmr_1d_spectra_options$solvent==solvent]]
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
      nmr_1d_spectra=nmr_1d_spectra[rownames(nmr_1d_spectra_options)[nmr_1d_spectra_options$ph<=ph[2]&nmr_1d_spectra_options$ph>=ph[1]]]
    }
    else if (length(ph)==1){
      nmr_1d_spectra=nmr_1d_spectra[rownames(nmr_1d_spectra_options)[nmr_1d_spectra_options$ph==ph]]
    }
  }
  
  #TEMPERATURE
  #25, 50, NA
  if (!is.null(temperature)){
    if(temperature%in%c(25,50)){
      nmr_1d_spectra=nmr_1d_spectra[rownames(nmr_1d_spectra_options)[nmr_1d_spectra_options$temperature==temperature]]
    }
    else{
      print("No temperature in library with that value.")
    }
  }
  
  return(nmr_1d_spectra[!is.na(names(nmr_1d_spectra))])
}





#' FIND THE BEST CORRELATION VALUE FOR THE CLUSTERING OF PEAKS
#' 
#' Takes a dataset and calculates the optimum correlation value, which leads to the maximum
#' number of clusters of the variables.
#' 
#' @param dataset_orig List representing the dataset from an nmr peaks metabolomics experiment.
#' @param CMETH Correlation method used to cluster the variables. Defaults to "pearson".
#' 
#' @return Value of the optimal correlation.
find_corr <- function(dataset_orig, CMETH='pearson') {
  
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
  CVAL_MAX<-0.999
  CVAL_STEP<-0.001
  CVAL_OPT_MIN<-0
  CVAL_OPT_MAX<-0
  CVAL_OPT<-0
  NBCLUSTERS_MAX<-0
  CVAL_CRIT<-0
  MAXSIZE<-40
  
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
#' 
#' Takes a dataset and performs clustering of variables, according to a correlation. The variables
#' will be separated into different clusters, according to a minimum correlation between variables.
#' Each cluster will correspond to a metabolite.
#' 
#' @param dataset_orig List representing the dataset from an nmr peaks metabolomics experiment.
#' @param CMETH Correlation method used to cluster the variables. Defaults to "pearson".
#' @param CVAL Minimum correlation between variables, so they can belong to the same cluster.
#' @param MVER Minimum number of variables in each cluster. Only the clusters with at least MVER variables will be returned.
#' 
#' @return List with the formed clusters
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
#' 
#' @param dataset List representing the dataset from an nmr peaks metabolomics experiment.
#' @param ppm.tol ppm tolerance when matching reference peaks to the dataset peaks.
#' @param clust.method Correlation method to use in the formation of clusters. Defaults to "pearson".
#' @param clust.treshold Minimum correlation between variables to form clusters. If not given, the function calculates the optimum value (value that leads to the greater number of clusters).
#' @param clust.peaks.min Minimum number of variables in each cluster. Only the clusters with at least clust.peaks.min variables will be considered.
#' @param freq Frequency of reference spectra. See documentation on \code{\link{choose_nmr_references}} function for further details.
#' @param nucl Atomic nuclei of reference spectra, either "1H" or "13C".
#' @param solv Solvent. See documentation on \code{\link{choose_nmr_references}} function for further details.
#' @param pH pH. See documentation on \code{\link{choose_nmr_references}} function for further details.
#' @param temp Temperature. See documentation on \code{\link{choose_nmr_references}} function for further details.
#' 
#' @return
#' List with of the results for each cluster. For each cluster:
#' \describe{
#'    \item{cluster.peaks}{The peaks of the cluster.}
#'    \item{summary}{Scores of each top 5 reference metabolite matched with the cluster.}
#'    \item{metabolites.matched}{List with full results for each top 5 reference metabolite matched.}
#' }
#' For each reference metabolite matched:
#' \describe{
#'    \item{score}{Jaccard index score.}
#'    \item{matched_peaks_ref}{Matched peaks of the reference spectra.}
#'    \item{matched_peaks_clust}{Matched peaks of the cluster.}
#'    \item{reference_peaks}{All the reference spectra peaks.}
#' }
#' 
#' @examples
#' 
#' 
#' @export
nmr_identification <- function(dataset, ppm.tol=0.03,
                               clust.method='pearson', clust.treshold=NULL, clust.peaks.min=2,
                               freq=500, nucl="1H", solv=NULL, pH=NULL, temp=NULL){
  
  #Choose reference metabolites:
  cat("Getting reference metabolites...\n")
  references=choose_nmr_references(frequency=freq, nucleus=nucl, solvent=solv, ph=pH,
                                   temperature=temp)
  
  #Get names of the reference metabolites:
  #split_ref=strsplit(names(references), split="_")
  #ref_names=c()
  #for (names in split_ref){
  #  ref_names=c(ref_names, names[1])
  #}
  data(nmr_1d_spectra_options, package="specmine")
  ref_names=as.character(nmr_1d_spectra_options[names(references), "met_id"])
  
  #Find best correlation value for the clusters correlation treshold, if a value is not given:
  if (is.null(clust.treshold)){
    cat("Getting best correlation value...\n")
    clust.treshold=find_corr(dataset, CMETH=clust.method)
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
    for (ref_id in 1:length(references)){
      res_ref=jaccard_index(clusters[[clust_id]], references[[ref_id]], PPMTOL=ppm.tol)
      
      score_clust=c(score_clust, res_ref$score)
      full_res[[ref_names[ref_id]]]=res_ref
    }
    #Store only the top 5 matches of the cluster and the cluster peaks
    names(score_clust)=ref_names
    score_clust=sort(score_clust, decreasing=T)[1:5]
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


