
metabolights_studies_list=function(){
  if (!requireNamespace("RCurl", quietly = TRUE)) {
    stop("Package RCurl needed for this function to work. Please install it: install.packages('RCurl')",
         call. = FALSE)
  }
  
  ftp_base="ftp://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/"
  listStudies_character=RCurl::getURL(ftp_base, dirlistonly=TRUE)
  listStudies_vector=strsplit(listStudies_character, "\n")[[1]]
  return(listStudies_vector)
}



get_files_list_per_assay=function(studyID){#, directory){
  ftp_base=paste("ftp://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/", studyID, "/", sep="")
  
  assays_files=list()
  
  i_file=readLines(paste(ftp_base, "i_Investigation.txt", sep=""))
  #curl::curl_download(paste(ftp_base, "i_Investigation.txt", sep=""), paste(directory, "i_Investigation.txt", sep="/"))
  #i_file=readLines(paste(directory, "i_Investigation.txt", sep="/"))
  
  assays_names_line=grep("^Study Assay File Name", i_file, value=TRUE)
  assays_names=grep("^a_", strsplit(assays_names_line, "\"")[[1]], value=TRUE)
  
  assay_files=paste(ftp_base, assays_names, sep="")
  #for (i in 1:length(assays_names)) curl::curl_download(paste(ftp_base, assays_names[i], sep=""), paste(directory, assays_names[i], sep="/"), quiet=F)
  #assay_files=paste(directory, assays_names, sep="/")
  
  
  for (i in 1:length(assay_files)){
    assay=read.table(assay_files[i], header=TRUE)
    if("Free.Induction.Decay.Data.File"%in%colnames(assay)) files_column="Free.Induction.Decay.Data.File"
    else files_column="Raw.Spectral.Data.File"
    
    assay_info=assay[,c("Sample.Name", files_column)]
    colnames(assay_info)=c("Samples", "Files")
    
    assays_files[[paste("Assay", i)]]=assay_info
  }
  
  return(assays_files)
}



get_metabolights_study_samples_files=function(studyID, assay, directory){
  assays_in_study=get_files_list_per_assay(studyID)#, directory)
  
  samples_files=assays_in_study[[assay]]
  colnames(samples_files)=NULL
  write.csv(samples_files, paste(directory, "samples_files.csv", sep="/"), row.names=F)
}



get_metabolights_study_files_assay=function(studyID, assay, directory){
  if (!requireNamespace("curl", quietly = TRUE)) {
    stop("Package curl needed for this function to work. Please install it: install.packages('curl')",
         call. = FALSE)
  }
  
  files_per_assay=get_files_list_per_assay(studyID)#, directory)
  files_to_download=as.character(files_per_assay[[assay]][,"Files"])
  
  ftp_base=paste("ftp://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/", studyID, "/", sep="")
  files_to_download_paths=paste(ftp_base, files_to_download, sep="")
  
  files_dest=paste(directory, assay, "data", files_to_download, sep="/")
  if (!dir.exists(paste(directory, assay, sep="/"))){
    dir.create(paste(directory, assay, sep="/"))
    dir.create(paste(directory, assay, "data", sep="/"))
  } 
  else if (!dir.exists(paste(directory, assay, "data", sep="/"))){
    dir.create(paste(directory, assay, "data", sep="/"))
  }
  
  for (i in 1:length(files_to_download)) curl::curl_download(files_to_download_paths[i], files_dest[i], quiet=FALSE)
}



get_metabolights_study_metadata_assay=function(studyID, assay, directory){
  
  #Get factor names:
  ftp_base=paste("ftp://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/", studyID, "/", sep="")
  
  i_file=readLines(paste(ftp_base, "i_Investigation.txt", sep=""))
  #curl::curl_download(paste(ftp_base, "i_Investigation.txt", sep=""), paste(directory, "i_Investigation.txt", sep="/"))
  #i_file=readLines(paste(directory, "i_Investigation.txt", sep="/"))
  
  factors_names=strsplit(grep("^Study Factor Name", i_file, value=TRUE), "\"")[[1]]
  factors_names_vec=grep("\t", factors_names, value=TRUE, invert=TRUE)
  
  factors=c()
  for (factor in factors_names_vec){
    factors=c(factors, paste("Factor.Value.", paste(strsplit(factor, " ")[[1]], collapse="."), ".", sep=""))
  }
  
  #Get samples in the assay:
  files_per_assay=get_files_list_per_assay(studyID)#, directory)
  samples_in_assay=as.character(files_per_assay[[assay]][,1])
  
  #Get file with the metadata information:
  sample_file_o=strsplit(grep("^Study File Name", i_file, value=TRUE), "\"")[[1]][2]
  sample_file=paste(ftp_base, sample_file_o, sep="")
  
  sample_metadata=read.table(sample_file, header=TRUE)
  #curl::curl_download(sample_file, paste(directory, sample_file_o, sep="/"))
  #sample_metadata=read.table(paste(directory, sample_file_o, sep="/"), header=T)
  
  metadata=as.data.frame(sample_metadata[,factors])
  rownames(metadata)=as.character(sample_metadata$Sample.Name)
  metadata=as.data.frame(metadata[samples_in_assay,])
  colnames(metadata)=factors_names_vec
  rownames(metadata)=as.character(samples_in_assay)
  
  metadata_filepath=paste(directory, "/metadata", assay, ".csv", sep="")
  write.csv(metadata, metadata_filepath)
}



get_metabolights_study=function(studyID, directory){
  if (!requireNamespace("curl", quietly = TRUE)) {
    stop("Package curl needed for this function to work. Please install it: install.packages('curl')",
         call. = FALSE)
  }
  
  message("Getting assays...\n")
  assays_in_study=get_files_list_per_assay(studyID)#, directory)
  cat("Done.\n")
  
  for (assay in 1:length(assays_in_study)){
    message("Assay ", assay)
    message("\nGetting files from assay ", assay, "...")
    get_metabolights_study_files_assay(studyID, assay, directory)
    message("\nDone.")
    message("\nGetting metadata from assay ", assay, "...")
    get_metabolights_study_metadata_assay(studyID, assay, directory)
    message("\nDone")
    message("\nGetting samples_files file from assay ", assay, "...")
    samples_files=assays_in_study[[assay]]
    colnames(samples_files)=NULL
    write.csv(samples_files, paste(directory, "samples_files.csv", sep="/"), row.names=FALSE)
    message("\nDone")
  }
}
