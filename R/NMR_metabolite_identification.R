library(specmine.datasets)


nmr_identification = function(dataset, ppm.tol, frequency_scores, solvent_scores, organism_scores,
                              method='Match_uniq', per.sample=FALSE, tresh_zero=0, alpha=10e-4){
  
  if(!method%in%c('Match_uniq', 'Hyper', 'Hyper_uniq')) stop('Invalid method. Valid methods: Match_uniq, Hyper or Hyper_uniq')
  
  if(is.null(frequency_scores)|is.null(solvent_scores)|is.null(organism_scores)){
    stop('Arguments missing: frequency_scores, solvent_score or organism_scores. Please check help page of the function to know how to set them')
  }
  
  res=list()
  
  #UNIQUENESS SCORES:
  if(method%in%c('Match_uniq','Hyper_uniq')){
    message("Calculating uniqueness scores...")
    uniq_scores=uniqueness_scores(0)
    message("Done\n")
  }
  else uniq_scores = NULL
  #Total number of peaks in library:
  if(method%in%c('Hyper','Hyper_uniq')) n_peaks_total = length(unique(unlist(spectra_list)))
  else n_peaks_total = NULL
  
  ###GETTING COMPOUNDS FOR ORGANISM SCORES
  message("Getting compounds for organism scores...\n")
  all_orgs=get_OrganismsCodes()
  cpds_groups=list()
  for(group in names(organism_scores)){
    if (!group%in%c("other", "not_in_kegg")) message("-- Getting compounds from ", group, "\n")
    if(is_organism(group)) cpds_groups[[group]]=compounds_in_organism(group)
    else if (is_group(group, all_orgs)) cpds_groups[[group]]=compounds_in_group(group, all_orgs)
  }
  message('Done\n')
  
  #PERFORM IDENTIFICATION:
  if(!per.sample){
    samples_peaks = as.numeric(rownames(dataset$data))
    return(identification_nmr_peaks(samples_peaks, method, ppm.tol, frequency_scores, solvent_scores,
                                    organism_scores, cpds_groups,
                                    uniq_scores, n_peaks_total, alpha))
  }
  res = list()
  for(samp in colnames(dataset$data)){
    message('-SAMPLE', samp, '\n')
    sample_peaks=as.numeric(names(dataset$data[dataset$data[,samp]>tresh_zero,1]))
    res[[samp]] = identification_nmr_peaks(sample_peaks, method, ppm.tol, frequency_scores, solvent_scores,
                                           organism_scores, cpds_groups,
                                           uniq_scores,n_peaks_total, alpha)
  }
  return(res)
}



identification_nmr_peaks = function(sample_peaks, method, ppm.tol, frequency_scores, solvent_scores,
                                    organism_scores, cpds_groups,
                                    uniq_scores=NULL, n_peaks_total=NULL, alpha=10e-4){
  
  message("Matching samples to reference peaks...")
  if(method=="Hyper") res = Hyper_method(sample_peaks, ppm.tol, n_peaks_total, alpha)
  else if(method=="Hyper_uniq") res = Hyper_uniq_method(sample_peaks, ppm.tol, uniq_scores, n_peaks_total, alpha)
  else res = Match_uniq_method(sample_peaks, ppm.tol, uniq_scores)
  message("Done...\n")
  
  message("Calculating Final scores...")
  score=c()
  score_frequency=c()
  score_solvents=c()
  score_organisms=c()
  for(i in 1:dim(res$results_table)[1]){
    reference_id=res$results_table$SPCMNS[i]
    
    score_fr=score_freq(reference_id, frequency_scores)
    score_solv=score_solvent(reference_id, solvent_scores)
    score_org=score_organism(res$results_table$SPCMNM[i], cpds_groups, organism_scores)
    score_final=(res$results_table$match_score[i]+score_fr+score_solv+score_org)/4
    
    score=c(score, score_final)
    score_frequency=c(score_frequency, score_fr)
    score_solvents=c(score_solvents, score_solv)
    score_organisms=c(score_organisms, score_org)
  }
  
  if(method=='Hyper'){
    res$results_table = cbind(res$results_table[,c('SPCMNM','Name','SPCMNS')],
                              score,
                              res$results_table[,c('match_score')],
                              score_frequency, score_solvents, score_organisms,
                              res$results_table[,c('n.peaks.matched','detailed_results_id')])
    colnames(res$results_table)[5] = 'match_score'
  }
  else if(method=='Hyper_uniq'){
    res$results_table = cbind(res$results_table[,c('SPCMNM','Name','SPCMNS')],
                              score,
                              res$results_table[,c('match_score', 'hypergeometric_score', 'uniqueness_score')],
                              score_frequency, score_solvents, score_organisms,
                              res$results_table[,c('n.peaks.matched','detailed_results_id')])
  }
  else{
    res$results_table = cbind(res$results_table[,c('SPCMNM','Name','SPCMNS')],
                              score,
                              res$results_table[,c('match_score', 'ratio', 'uniqueness_score')],
                              score_frequency, score_solvents, score_organisms,
                              res$results_table[,c('n.peaks.matched','detailed_results_id')])
  }
  colnames(res$results_table)[4] = 'Final_Score'
  res$results_table=res$results_table[order(res$results_table$Final_Score, decreasing=TRUE),]
  message("Done.\n")
  return(res)
}


Hyper_method = function(peaks, ppm.tol, n_peaks_total, alpha=10e-4){
  res = list()
  Name=c()
  spcmnm=c()
  ratio=c()
  uniqueness=c()
  match_score=c()
  p.fdr=c()
  n.peaks.matched=c()
  
  spcmns=c()
  j=0
  for (reference_id in names(spectra_list)){
    reference=spectra_list[[reference_id]]
    
    res.ref=match_Hyper(peaks, unique(reference$chemical.shift), n_peaks_total, ppm.tol)
    
    if(res.ref$matched_peaks!=0){
      spcmnm_ids=get_spcmnm_from_spcmns(reference_id)
      for (i in 1:length(spcmnm_ids)){
        spcmns=c(spcmns, reference_id)
        j=j+1
        id=paste(reference_id, j, sep="_")
        res$more_results[[id]]=res.ref[c("matched_peaks_ref", "matched_peaks_samp", "reference_peaks")]
        spcmnm=c(spcmnm, spcmnm_ids[i])
        Name=c(Name, codes[codes$SPCMNM==spcmnm_ids[i], "NAME"])
        match_score=c(match_score, res.ref$score)
        n.peaks.matched=c(n.peaks.matched, res.ref$matched_peaks)
      }
    }
  }
  
  p.fdr=p.adjust(match_score, "fdr")
  hypergeometric_score = rep(0, length(p.fdr))
  hypergeometric_score[p.fdr>alpha] = 0
  hypergeometric_score[p.fdr<alpha] = (1-(p.fdr[p.fdr<alpha]/alpha))
  
  res$results_table=data.frame(SPCMNM=spcmnm, Name, SPCMNS=spcmns, match_score=hypergeometric_score, n.peaks.matched, detailed_results_id=names(res$more_results), stringsAsFactors=FALSE)
  maintain_non_zeros = res$results_table$match_score!=0
  maintain_non_zeros_detailed = res$results_table$detailed_results_id[maintain_non_zeros]
  res$results_table = res$results_table[maintain_non_zeros,]
  res$more_results = res$more_results[maintain_non_zeros_detailed]
  
  
  return(res)
}



Hyper_uniq_method = function(peaks, ppm.tol, uniq_scores, n_peaks_total, alpha=10e-4){
  res = list()
  Name=c()
  spcmnm=c()
  ratio=c()
  uniqueness=c()
  match_score=c()
  p.fdr=c()
  n.peaks.matched=c()
  
  spcmns=c()
  j=0
  for (reference_id in names(spectra_list)){
    reference=spectra_list[[reference_id]]
    
    res.ref=match_Hyper(peaks, unique(reference$chemical.shift), n_peaks_total, ppm.tol)
    
    if(res.ref$matched_peaks!=0){
      spcmnm_ids=get_spcmnm_from_spcmns(reference_id)
      for (i in 1:length(spcmnm_ids)){
        spcmns=c(spcmns, reference_id)
        j=j+1
        id=paste(reference_id, j, sep="_")
        res$more_results[[id]]=res.ref[c("matched_peaks_ref", "matched_peaks_samp", "reference_peaks")]
        spcmnm=c(spcmnm, spcmnm_ids[i])
        Name=c(Name, codes[codes$SPCMNM==spcmnm_ids[i], "NAME"])
        match_score=c(match_score, res.ref$score)
        n.peaks.matched=c(n.peaks.matched, res.ref$matched_peaks)
      }
    }
  }
  
  p.fdr=p.adjust(match_score, "fdr")
  hypergeometric_score = rep(0, length(p.fdr))
  hypergeometric_score[p.fdr>alpha] = 0
  hypergeometric_score[p.fdr<alpha] = (1-(p.fdr[p.fdr<alpha]/alpha))
  
  spcmns_uniq_scores = c(NA, length(hypergeometric_score))
  Hyper_uniq_score = rep(0, length(hypergeometric_score))
  for(i in 1:length(hypergeometric_score)){
    if(hypergeometric_score[i]!=0){
      Hyper_uniq_score[i] = (hypergeometric_score[i] + uniq_scores[spcmns[i]]) / 2
      spcmns_uniq_scores[i] = uniq_scores[spcmns[i]]
    }
  }
  
  res$results_table=data.frame(SPCMNM=spcmnm, Name, SPCMNS=spcmns, match_score=Hyper_uniq_score, hypergeometric_score, uniqueness_score=spcmns_uniq_scores, n.peaks.matched, detailed_results_id=names(res$more_results), stringsAsFactors=FALSE)
  maintain_non_zeros = res$results_table$match_score!=0
  maintain_non_zeros_detailed = res$results_table$detailed_results_id[maintain_non_zeros]
  res$results_table = res$results_table[maintain_non_zeros,]
  res$more_results = res$more_results[maintain_non_zeros_detailed]
  
  return(res)
}



Match_uniq_method = function(peaks, ppm.tol, uniq_scores){
  res = list()
  Name=c()
  spcmnm=c()
  ratio=c()
  uniqueness=c()
  match_score=c()
  p.fdr=c()
  n.peaks.matched=c()
  
  spcmns=c()
  j=0
  for (reference_id in names(spectra_list)){
    reference=spectra_list[[reference_id]]
    
    res.ref=match_Match_uniq(peaks, unique(reference$chemical.shift), uniq_scores[[reference_id]], ppm.tol)
    
    if(res.ref$matched_peaks!=0){
      spcmnm_ids=get_spcmnm_from_spcmns(reference_id)
      for (i in 1:length(spcmnm_ids)){
        spcmns=c(spcmns, reference_id)
        j=j+1
        id=paste(reference_id, j, sep="_")
        res$more_results[[id]]=res.ref[c("matched_peaks_ref", "matched_peaks_samp", "reference_peaks")]
        spcmnm=c(spcmnm, spcmnm_ids[i])
        Name=c(Name, codes[codes$SPCMNM==spcmnm_ids[i], "NAME"])
        ratio=c(ratio, res.ref$ratio)
        uniqueness=c(uniqueness, res.ref$uniqueness)
        match_score=c(match_score, res.ref$score)
        n.peaks.matched=c(n.peaks.matched, res.ref$matched_peaks)
      }
    }
  }
  
  res$results_table=data.frame(SPCMNM=spcmnm, Name, SPCMNS=spcmns, match_score, ratio, uniqueness_score=uniqueness, n.peaks.matched, detailed_results_id=names(res$more_results), stringsAsFactors=FALSE)
  
  return(res)
}



match_Hyper = function(sample_peaks, refMetab_peaks, n_peaks_total, PPMTOL=0.03){
  l.samp=length(sample_peaks)
  l.ref=length(refMetab_peaks)
  
  res=list()
  
  #Save the peaks from the cluster and references that matched:
  matched_peaks_samp=c()
  matched_peaks_ref=c()
  
  #Sort sample and reference peaks:
  sample_peaks=sort(sample_peaks)
  refMetab_peaks=sort(refMetab_peaks)
  
  #Variables to calculate hypergeometric test:
  matched_peaks=0
  n_ref_peaks=l.ref
  n_notRef_peaks=n_peaks_total-n_ref_peaks
  n_samp_peaks=l.samp
  
  if (l.samp<l.ref){
    #For each sample peak, compare it to the reference peaks:
    last_pR=1
    for(pS in 1:length(sample_peaks)){
      if (last_pR>l.ref) break
      for(pR in last_pR:l.ref){
        inf_limit=refMetab_peaks[pR]-PPMTOL
        if(round(sample_peaks[pS],3) < round(inf_limit,3)) break
        sup_limit=refMetab_peaks[pR]+PPMTOL
        if(round(sample_peaks[pS],3) <= round(sup_limit,3)){
          #if(multiplets){
          #  n_multiplets = sum(refMetab_peaks == refMetab_peaks[pR])
          #  matched_peaks = matched_peaks + n_multiplets
          #}
          #else{
          matched_peaks = matched_peaks + 1
          last_pR=pR+1
          #}
          
          matched_peaks_samp=c(matched_peaks_samp, round(sample_peaks[pS],3))
          matched_peaks_ref=c(matched_peaks_ref, round(refMetab_peaks[pR],3))
          break
        }
      }
    }
  }
  else{
    #For each reference peak, compare it to the sample peaks:
    last_pS=1
    #if(multiplets) refids = which(!duplicated(refMetab_peaks))
    #else refids = 1:refMetab_peaks
    for(pR in 1:length(refMetab_peaks)){
      if (last_pS>l.samp) break
      for(pS in last_pS:l.samp){
        inf_limit=sample_peaks[pS]-PPMTOL
        if(round(refMetab_peaks[pR],3) < round(inf_limit,3)) break
        sup_limit=sample_peaks[pS]+PPMTOL
        if(round(refMetab_peaks[pR],3) <= round(sup_limit,3)){
          #if(multiplets) matched_peaks = matched_peaks + sum(refMetab_peaks == refMetab_peaks[pR])
          #else
          matched_peaks = matched_peaks + 1
          last_pS=pS+1
          
          matched_peaks_samp=c(matched_peaks_samp, round(sample_peaks[pS],3))
          matched_peaks_ref=c(matched_peaks_ref, round(refMetab_peaks[pR],3))
          break
        }
      }
    }
  }
  
  #Calculate hypergeometric test:
  ht=dhyper(matched_peaks, n_ref_peaks, n_notRef_peaks, n_samp_peaks)
  
  #Store results:
  res$matched_peaks_ref=matched_peaks_ref
  res$matched_peaks_samp=matched_peaks_samp
  res$reference_peaks=refMetab_peaks
  res$score=ht
  res$matched_peaks=matched_peaks
  
  return(res)
}



match_Match_uniq = function(sample_peaks, refMetab_peaks, uniq_score, PPMTOL=0.03){
  l.samp=length(sample_peaks)
  l.ref=length(refMetab_peaks)
  
  res=list()
  
  #Save the peaks from the cluster and references that matched:
  matched_peaks_samp=c()
  matched_peaks_ref=c()
  
  #Sort sample and reference peaks:
  sample_peaks=sort(sample_peaks)
  refMetab_peaks=sort(refMetab_peaks)
  
  matched_peaks=0
  n_ref_peaks=l.ref
  n_samp_peaks=l.samp
  if (l.samp<l.ref){
    #For each sample peak, compare it to the reference peaks:
    last_pR=1
    for(pS in 1:length(sample_peaks)){
      if (last_pR>l.ref) break
      for(pR in last_pR:l.ref){
        inf_limit=refMetab_peaks[pR]-PPMTOL
        if(round(sample_peaks[pS],3) < round(inf_limit,3)) break
        sup_limit=refMetab_peaks[pR]+PPMTOL
        if(round(sample_peaks[pS],3) <= round(sup_limit,3)){
          #if(multiplets){
          #  n_multiplets = sum(refMetab_peaks == refMetab_peaks[pR])
          #  matched_peaks = matched_peaks + n_multiplets
          #}
          #else{
          matched_peaks = matched_peaks + 1
          last_pR=pR+1
          #}
          
          matched_peaks_samp=c(matched_peaks_samp, round(sample_peaks[pS],3))
          matched_peaks_ref=c(matched_peaks_ref, round(refMetab_peaks[pR],3))
          break
        }
      }
    }
  }
  else{
    #For each reference peak, compare it to the sample peaks:
    last_pS=1
    #if(multiplets) refids = which(!duplicated(refMetab_peaks))
    #else refids = 1:refMetab_peaks
    for(pR in 1:length(refMetab_peaks)){
      if (last_pS>l.samp) break
      for(pS in last_pS:l.samp){
        inf_limit=sample_peaks[pS]-PPMTOL
        if(round(refMetab_peaks[pR],3) < round(inf_limit,3)) break
        sup_limit=sample_peaks[pS]+PPMTOL
        if(round(refMetab_peaks[pR],3) <= round(sup_limit,3)){
          #if(multiplets) matched_peaks = matched_peaks + sum(refMetab_peaks == refMetab_peaks[pR])
          #else
          matched_peaks = matched_peaks + 1
          last_pS=pS+1
          
          matched_peaks_samp=c(matched_peaks_samp, round(sample_peaks[pS],3))
          matched_peaks_ref=c(matched_peaks_ref, round(refMetab_peaks[pR],3))
          break
        }
      }
    }
  }
  
  #Calculate Match_uniq score:
  if(matched_peaks==0) score = 0
  else score = ((matched_peaks/n_ref_peaks) + uniq_score) / 2
  
  #Store results:
  res$matched_peaks_ref=matched_peaks_ref
  res$matched_peaks_samp=matched_peaks_samp
  res$reference_peaks=refMetab_peaks
  res$score=score
  res$ratio=matched_peaks/n_ref_peaks
  if(matched_peaks!=0) res$uniqueness=uniq_score
  else res$uniqueness=NA
  res$matched_peaks=matched_peaks
  
  return(res)
}


uniqueness_scores = function(ppm_tolerance){
  
  all_peaks = c()
  for(spec in names(spectra_list)) all_peaks = c(all_peaks, unique(round(spectra_list[[spec]]$chemical.shift, digits=3)))
  peaks_table = table(all_peaks)
  
  uniq_values_peaks = c()
  for(peak in unique(all_peaks)){
    overlap = sum(peaks_table[as.character(all_peaks[all_peaks>=peak-ppm_tolerance & all_peaks<=peak+ppm_tolerance])])
    uniq_values_peaks = c(uniq_values_peaks, 1/overlap)
  }
  names(uniq_values_peaks) = as.character(unique(all_peaks))
  
  res = c()
  for(spec in names(spectra_list)){
    peaks = as.character(unique(round(spectra_list[[spec]]$chemical.shift, digits=3)))
    uniq_score_spec = mean(uniq_values_peaks[peaks])
    res = c(res, uniq_score_spec)
  }
  names(res) = names(spectra_list)
  
  return(res)
}




###DATA CONVERSION
##SPCMNS-->SPCMNM
get_spcmnm_from_spcmns=function(spcmns){
  metabolites=unlist(strsplit(spectra_options[spectra_options$SPCMNS==spcmns, "SPCMNM"], "; "))
  return(metabolites)
}

##SPCMNM-->KEGG
convert_spcmnm_to_kegg=function(spcmnm){
  keggs=unlist(strsplit(codes$KEGG[codes$SPCMNM==spcmnm], "; ")[[1]])
  return(keggs)
}

##CHEBI-->SPCMNM
convert_chebi_to_spcmnm=function(chebi){
  spcmnms=c()
  for(i in 1:length(codes$CHEBI)){
    chebi_ref=codes$CHEBI[i]
    if (chebi %in% unlist(strsplit(chebi_ref, "; ")[[1]])) spcmnms=c(spcmnms, codes$SPCMNM[i])
  }
  return(spcmnms)
}



###SCORES
##Frequency
score_freq=function(spcmns, scores=list('400'=0.5, '500'=1, '600'=0.5, '700'=0.5)){
  score=scores[[as.character(spectra_options[spectra_options$SPCMNS==spcmns, "FREQUENCY"])]]
  return(score)
}

##Solvent
score_solvent=function(spcmns, scores=list(CD3OD=1, D2O=0.8, Water=0.8, CDCl3=0.6, 'Acetone-d6'=0.6, 'Acetone'=0.6,
                                           "DMSO-d6"=0.4, "100%_DMSO"=0.4, "5%_DMSO"=0.4, C=0.2, C6D6=0.2,
                                           CD3CN=0.2, C2D2Cl4=0.2, CD2Cl2=0.2, CDC3OD=0.2, Ethanol=1)){
  score=scores[[spectra_options[spectra_options$SPCMNS==spcmns, "SOLVENT"]]]
  return(score)
}

##organism
score_organism=function(spcmnm, cpds_groups, scores=list('hsa'=1, 'Mammals'=0.8, 'Vertebrates'=0.7, 'Animals'=0.6, 'Eukaryotes'=0.3, 'other'=0.1, 'not_in_kegg'=0.6)){
  #To specify organism: Kegg code
  #Availale groups: "Eukaryotes", "Animals", "Vertebrates", "Mammals", "Birds", "Reptiles", "Amphibians", "Fishes", "Arthropods", "Insects",
  #"Nematodes", "Mollusks", "Cnidarians", "Plants", "Eudicots", "Monocots", "Green algae", "Red algae", "Fungi", "Protist", "Prokaryotes",
  #"Bacteria", "Archaea"
  
  cpds=convert_spcmnm_to_kegg(spcmnm)
  
  if(is.null(cpds)) return(scores[['not_in_kegg']])
  else{
    possible_scores=c(scores[['other']])
    for(cpd in cpds){
      for(group in names(scores)){
        if(!group%in%c('other', 'not_in_kegg')) if(cpd%in%cpds_groups[[group]]) possible_scores=c(possible_scores, scores[[group]])
      }
    }
    score=max(possible_scores)
  }
  return(score)
}
