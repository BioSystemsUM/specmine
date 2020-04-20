

#################################################
#####READ NMR SPECTRA FILES IN BRUKER FORMAT#####
#################################################
##One Bruker directory corresponds to one bruker spectrum:
###spectrum_directory/
#####-pdata/
#######-1/
#########-1r
#########-procs
#########-other files
#####-other files
readBruker<-function(BrukerDataDir)
{
  ## written by Dr. Jie Hao, Imperial College London
  ##Taken from nmrML github repository
  warnDef<-options("warn")$warn
  warnRead<-options(warn = -1)
  datapath<-BrukerDataDir
  
  ## read in bruker spectra
  ## find the data files
  ppm <- NULL
  pfile <-list.files(path = datapath, pattern = "^procs$", all.files = FALSE,full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  rfile <-list.files(path = datapath, pattern = "^1r$", all.files = FALSE,full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  L<-length(pfile)
  Lr<-length(rfile)
  sa <- NULL
  snam <- NULL
  if (L==0 || Lr==0 || L!=Lr)
  {
    return (cat("Bruker file does not exist in datapath, or other problems with bruker files...\n"))
  } else {
    for (i in 1:L)
    {    
      con  <- file(pfile[i], open = "r")
      aLine <- readLines(con, n = -1, warn = FALSE)
      myV <- strsplit(aLine, "=")    
      close(con)
      
      for (j in 1:length(myV))
      {
        if (match("##$OFFSET",myV[[j]][1],nomatch = 0))
        {    offset <- as.numeric(myV[[j]][2]);
        }
        if (match("##$SW_p",myV[[j]][1],nomatch = 0))
        {    sw <- as.numeric(myV[[j]][2]);
        }
        if (match("##$SF",myV[[j]][1],nomatch = 0))
        {
          sf <- as.numeric(myV[[j]][2]);
        }
        if (match("##$SI",myV[[j]][1],nomatch = 0))
        {  
          si <- as.numeric(myV[[j]][2]);
        }
        if (match("##$BYTORDP",myV[[j]][1],nomatch = 0))
        {    bytordp <- as.numeric(myV[[j]][2]);
        }
        if (match("##$NC_proc",myV[[j]][1],nomatch = 0))
        {
          ncproc <- as.numeric(myV[[j]][2]);
        }
      }
      
      if (bytordp==0){machine_format =  "little"}
      else {machine_format = "big"}
      
      s<-readBin(rfile[i], what="int",70000, size = 4, signed = T, endian =machine_format)
      s<- ((2^ncproc)* s)
      nspec <- length(s)
      
      tmpppm <- ppm
      
      swp <- sw/sf
      dppm <- swp/(nspec-1)
      ppm<-offset
      ppm<-seq(offset,(offset-swp),by=-dppm)
      
      ## interpolation
      if (!is.null(tmpppm))
      {
        if (length(tmpppm) != length(ppm))
        {
          sinter <- approx(ppm, s, xout = tmpppm)
          s <- sinter$y
          s[is.na(s)]<-0
          ppm <- tmpppm
        }
      }
      
      sa<- cbind(sa,s)
      ## find corresponding title
      stitle<-paste(substr(rfile[i],1,nchar(rfile[i])-2),"title",sep="")
      if (!file.exists(stitle))
        stitle<-paste(substr(rfile[i],1,nchar(rfile[i])-2),"TITLE",sep="")
      if (file.exists(stitle))
      {
        if (!file.info(stitle)$size == 0)
        {
          con<-file(stitle,open="r")
          ntem <- readLines(con, n = 1, warn = FALSE)
          close(con)
        } else {
          sT <- strsplit(rfile[i], "/")
          sTitle <-sT[[1]]         
          lsT<-length(sTitle)
          if (lsT>4)
            ntem<-paste(sTitle[lsT-4],"_",sTitle[lsT-3],"_",sTitle[lsT-1],sep="")
          else if (lsT>3)
            ntem<-paste(sTitle[lsT-3],"_",sTitle[lsT-1],sep="")
          else if (lsT>=1)
            ntem<-paste(sTitle[lsT-1],sep="")
          else
            ntem<-i
        }
      } else {
        sT <- strsplit(rfile[i], "/")
        sTitle <-sT[[1]]         
        lsT<-length(sTitle)
        if (lsT>4)
          ntem<-paste(sTitle[lsT-4],"_",sTitle[lsT-3],"_",sTitle[lsT-1],sep="")
        else if (lsT>3)
          ntem<-paste(sTitle[lsT-3],"_",sTitle[lsT-1],sep="")
        else if (lsT>=1)
          ntem<-paste(sTitle[lsT-1],sep="")
        else
          ntem<-i
      }
      snam<- cbind(snam, ntem)            
    }
  }
  snam <- cbind("ppm", snam)
  sa <- cbind(ppm,sa)
  colnames(sa)<- snam
  warnRead<-options(warn = warnDef)
  return (sa)
}



read_Bruker_files=function(bruker_directory, metadata_file=NULL,
                           m.header_col=T, m.header_row=T, m.sep=",",
                           samples.names=NULL, zipped=T,
                           description="", label.x="ppm", label.values="intensity"){
  
  #Get data paths of the files to read:
  files_study=list.files(bruker_directory, include.dirs=T, full.names=T, recursive=T)
  if (zipped){
    folders_zip_data=grep("[.]zip$", files_study, value=T)
    if(length(folders_zip_data)==0) stop("No zip files in directory ", bruker_directory)
    #Read each spectra files:
    temp_directory="~/temp" #tempdir()
    dir.create(temp_directory)
    for (dir in folders_zip_data){
      x=tail(unlist(strsplit(dir, "/")), n=1)
      folder_name=unlist(strsplit(x, "[.]zip"))
      x=unzip(dir, exdir=paste(temp_directory, folder_name, sep="/"))
    }
    files_bruker=list.files(temp_directory, include.dirs=T, full.names=T, recursive=T)
  }
  else files_bruker=files_study 
  directories_to_read=grep("/pdata/1$", files_bruker, value=T)
  
  #Get samples names:
  samples_dirs_list=c() #Name of the directories of each sample
  samples_names_list=c() #Name of each sample
  #Sample name in index 1 of samples_names_list corresponds to sample directory in index 1 of samples_dirs_list ^
  
  if (!is.null(samples.names)){
    samples_names_df=read.csv(samples.names, header=F, stringsAsFactors=F)
    if(length(grep("zip$", samples_names_df[[2]], value=T))!=0){
      samples_names_df[[2]]=substr(samples_names_df[[2]], 1, nchar(as.character(samples_names_df[[2]]))-4)
    }
    samples_dirs_list=samples_names_df[,2]
    samples_names_list=samples_names_df[,1]
  }
  else{
    for(dir in directories_to_read){
      n=tail(strsplit(dir, "/")[[1]], n=4)
      n=n[1]
      samples_names_list=c(samples_names_list, n)
      samples_dirs_list=c(samples_dirs_list, n)
    }
  }
  if (length(directories_to_read) != length(samples_names_list)) stop("File structure is incorrect. Maybe you have more than one spectrum in a sample folder or the samples folders inside each .zip have repeated names. Please look at the documentation.")
  
  #Read metadata file, if given:
  if(!is.null(metadata_file)){
    cat("Reading Metadata file\n")
    metadata=specmine::read_metadata(metadata_file, header.col=m.header_col, header.row=m.header_row, sep=m.sep)
  }
  else metadata=NULL
  
  
  spectra.list=list()
  #Read each spectra files:
  data_matrix=NULL
  for (i in 1:length(directories_to_read)){
    dir=directories_to_read[i]
    sample_name=samples_names_list[match(tail(unlist(strsplit(dir, "/")), n=4)[1], samples_dirs_list)]
    
    cat("Reading sample ", sample_name, " in ", dir)
    cat("\n")
    
    matrix_sample=readBruker(dir)
    if (!is.null(matrix_sample)){
      spectra.list[[sample_name]]=data.frame(ppm=matrix_sample[,1], intensity=matrix_sample[,2])
      j=order(spectra.list[[sample_name]][,"ppm"])
      spectra.list[[sample_name]]=spectra.list[[sample_name]][j,]
      spectra.list[[sample_name]][,"ppm"]=round(spectra.list[[sample_name]][,"ppm"], 3)
    }
  }
  
  cat("Creating dataset (this may take a while)\n")
  dataset=specmine::dataset_from_peaks(spectra.list, type="nmr-spectra", metadata=metadata, description=description)
  
  if(zipped) unlink(temp_directory, recursive=T)
  
  cat("Done.")
  
  return(dataset)
}






#################################################
#####READ NMR SPECTRA FILES IN VARIAN FORMAT#####
#################################################

is_windows <- function() {
  identical(.Platform$OS.type, "windows")
}

read_varian_spectrum_raw=function(directory, zero_filling=T, apodization=T){
  
  files=list.files(directory, recursive=T, full.names = T)
  fid_file=grep(".*fid$", files, value=T)
  procpar_file=grep(".*procpar.*$", files, value=T)
  
  if(!reticulate::py_available()){
    p=reticulate::py_discover_config()
    python_versions=p$python_versions
    python_path=""
    if(is_windows() & length(python_versions)>0) python_path=python_versions[length(python_versions)]
    else python_path=grep(".*python3.*", p$python_versions, value=T)
    if(nchar(python_path)==0) stop("Python version 3 must be installed")
    reticulate::use_python(python_path)
    reticulate::py_config()
  }
  
  varian_raw_py=system.file("read_varian_spec_raw.py", package="specmine")
  env=new.env()
  reticulate::source_python(varian_raw_py,  envir = env)
  
  spec=env$read_varian_spec_raw(directory, fid_file, procpar_file, zero_filling=zero_filling, apodization=apodization)
  names(spec)=c("ppm", "intensity")
  return(spec)
}




read_varian_spectra_raw=function(varian_spectra_directory,
                                 metadata_file=NULL, m.header_col=T, m.header_row=T, m.sep=",",
                                 samples.names=NULL, zero_filling=T, apodization=T, zipped=T,
                                 description="", label.x="ppm", label.values="intensity"){
  
  spectra.list=list()
  
  #Get data paths of the files to read:
  if (zipped){
    files_study=list.files(varian_spectra_directory, include.dirs=T, full.names=T)
    folders_zip_data=grep("fid.zip$", files_study, value=T)
    #Read each spectra files:
    temp_directory="~/temp" #tempdir()
    dir.create(temp_directory)
    for (dir in folders_zip_data){
      unzip(dir, exdir=temp_directory)
    }
    directories_to_read=grep(".*[.]fid", list.dirs(temp_directory, full.names=T), value=T)
  }
  else directories_to_read=grep(".*[.]fid", list.dirs(varian_spectra_directory, full.names=T), value=T)
  
  #Get samples names:
  if (!is.null(samples.names)){
    samples_names=read.csv(samples.names, header=F)
    if(length(grep("zip$", samples_names[[2]], value=T))!=0){
      samples_names[[2]]=substr(samples_names[[2]], 1, nchar(as.character(samples_names[[2]]))-4)
    }
  }
  else{
    samples_names=c()
    for(dir in directories_to_read){
      n=tail(strsplit(dir, "/")[[1]], n=1)
      n=strsplit(n, ".fid")[[1]]
      samples_names=c(samples_names, n)
    }
  }
  
  #Get spectra, processed accordingly:
  cat("Reading files:\n")
  data_matrix_column_names=c()
  for (i in 1:length(directories_to_read)){
    dir=directories_to_read[i]
    cat(dir)
    cat("\n")
    if (!is.null(samples.names)) sample_name=samples_names[1][samples_names[2]==tail(strsplit(dir, "/")[[1]], n=1)]
    else sample_name=samples_names[i]
    spectrum=read_varian_spectrum_raw(dir, zero_filling=zero_filling, apodization=apodization)
    spectra.list[[sample_name]]=data.frame(spectrum)
    data_matrix_column_names=c(data_matrix_column_names, sample_name)
  }
  
  data_matrix=data.frame(spectra.list[[1]]$intensity, row.names=as.character(spectra.list[[1]]$ppm))
  for (i in 2:length(spectra.list)){
    data_matrix=cbind(data_matrix, spectra.list[[i]]$intensity)
  }
  colnames(data_matrix)=data_matrix_column_names
  
  if(!is.null(metadata_file)) metadata=specmine::read_metadata(metadata_file, header.col=m.header_col, header.row=m.header_row, sep=m.sep)
  else metadata=NULL
  
  dataset=specmine::create_dataset(data_matrix, type="nmr-spectra", metadata=metadata, description=description, label.x=label.x, label.values=label.values)
  
  if(zipped) unlink(temp_directory, recursive=T)
  
  return(dataset)
  
}






##########################
#####DETECT NMR PEAKS#####
##########################

detect_nmr_peaks=function(spectrum.list, baseline_treshold=50000){
  
  res=list()
  
  x=t(as.matrix(spectrum.list$intensity))
  peaks_detected_idx=speaq::detectSpecPeaks(x, verbose=F, baselineThresh=baseline_treshold)
  
  res$ppm=as.numeric(spectrum.list$ppm[peaks_detected_idx[[1]]])
  res$intensity=as.numeric(spectrum.list$intensity[peaks_detected_idx[[1]]])
  rownames(res)=NULL
  
  return(res)
}



detect_nmr_peaks_from_dataset=function(dataset, baseline_tresh=50000,
                                       ap.method="own", ap.samp.classes=1, ap.step=0.03){
  
  detected_peaks_list=list()
  for (spectrum in colnames(dataset$data)){
    
    spec=list()
    spec[["ppm"]]=rownames(dataset$data)
    spec[["intensity"]]=dataset$data[,spectrum]
    
    detected_peaks=detect_nmr_peaks(spec, baseline_treshold=baseline_tresh)
    
    cat("Spectrum ", spectrum, " has ", length(detected_peaks$ppm), " peaks.")
    cat("\n")
    
    detected_peaks_list[[spectrum]]=as.data.frame(detected_peaks)
    detected_peaks_list[[spectrum]]$ppm=round(detected_peaks_list[[spectrum]]$ppm, 2)
    
  }
  
  final_dataset=specmine::group_peaks(detected_peaks_list, type="nmr-peaks", metadata=dataset$metadata,
                                      method=ap.method, samp.classes=ap.samp.classes, step=ap.step,
                                      label.x=dataset$labels$x, label.values=dataset$labels$val, description=dataset$description)
  
  return(final_dataset)
  
}
