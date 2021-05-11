

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
  suppressWarnings({
    ## written by Dr. Jie Hao, Imperial College London
    ##Taken from nmrML github repository
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
      stop("Bruker file does not exist in datapath, or other problems with bruker files...\n")
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
        
        s<-readBin(rfile[i], what="int",70000, size = 4, signed = TRUE, endian =machine_format)
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
  })
  return (sa)
}



read_Bruker_files=function(bruker_directory, metadata_file=NULL,
                           m.header_col=TRUE, m.header_row=TRUE, m.sep=",",
                           samples.names=NULL, zipped=TRUE,
                           description="", label.x="ppm", label.values="intensity"){
  
  #Get data paths of the files to read:
  files_study=list.files(bruker_directory, include.dirs=T, full.names=TRUE, recursive=TRUE)
  if (zipped){
    folders_zip_data=grep("[.]zip$", files_study, value=TRUE)
    if(length(folders_zip_data)==0) stop("No zip files in directory ", bruker_directory)
    #Read each spectra files:
    temp_directory="~/temp" #tempdir()
    dir.create(temp_directory)
    for (dir in folders_zip_data){
      x=tail(unlist(strsplit(dir, "/")), n=1)
      folder_name=unlist(strsplit(x, "[.]zip"))
      x=unzip(dir, exdir=paste(temp_directory, folder_name, sep="/"))
    }
    files_bruker=list.files(temp_directory, include.dirs=TRUE, full.names=TRUE, recursive=TRUE)
  }
  else files_bruker=files_study 
  directories_to_read=grep("/pdata/1$", files_bruker, value=TRUE)
  
  #Get samples names:
  samples_dirs_list=c() #Name of the directories of each sample
  samples_names_list=c() #Name of each sample
  #Sample name in index 1 of samples_names_list corresponds to sample directory in index 1 of samples_dirs_list ^
  
  if (!is.null(samples.names)){
    samples_names_df=read.csv(samples.names, header=FALSE, stringsAsFactors=FALSE)
    if(length(grep("zip$", samples_names_df[[2]], value=TRUE))!=0){
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
    message("Reading Metadata file\n")
    metadata=specmine::read_metadata(metadata_file, header.col=m.header_col, header.row=m.header_row, sep=m.sep)
  }
  else metadata=NULL
  
  
  spectra.list=list()
  #Read each spectra files:
  data_matrix=NULL
  for (i in 1:length(directories_to_read)){
    dir=directories_to_read[i]
    sample_name=samples_names_list[match(tail(unlist(strsplit(dir, "/")), n=4)[1], samples_dirs_list)]
    
    message("Reading sample ", sample_name, " in ", dir)
    message("\n")
    
    matrix_sample=readBruker(dir)
    if (!is.null(matrix_sample)){
      spectra.list[[sample_name]]=data.frame(ppm=matrix_sample[,1], intensity=matrix_sample[,2])
      j=order(spectra.list[[sample_name]][,"ppm"])
      spectra.list[[sample_name]]=spectra.list[[sample_name]][j,]
      spectra.list[[sample_name]][,"ppm"]=round(spectra.list[[sample_name]][,"ppm"], 3)
    }
  }
  
  message("Creating dataset (this may take a while)\n")
  dataset=specmine::dataset_from_peaks(spectra.list, type="nmr-spectra", metadata=metadata, description=description)
  
  if(zipped) unlink(temp_directory, recursive=TRUE)
  
  message("Done.")
  
  return(dataset)
}

#################################################
#####READ 2D NMR SPECTRA FILES IN BRUKER FORMAT#####
#################################################
##One Bruker directory corresponds to one 2D bruker spectrum:
###spectrum_directory/
#####-pdata/
#######-1/
#########-2rr
#########-procs
#########-proc2s
#########-other files
#####-acqus
#####-acqu2s
#####-other files

#Uses function readBruker from mrbin package to read 2D Bruker files
read_Bruker_files_2d <- function(bruker_directory, metadata_file=NULL,
                                 m.header_col=T, m.header_row=T, m.sep=",",
                                 samples.names=NULL, zipped=T,
                                 description="", label.x="ppm", label.y = "ppm", label.values="intensity"){
  
  #Get data paths of the files to read:
  files_study <- list.files(bruker_directory, include.dirs=T, full.names=T)
  if (zipped){
    folders_zip_data <- grep("[.]zip$", files_study, value=T)
    if(length(folders_zip_data)==0) stop("No zip files in directory ", bruker_directory)
    #Read each spectra files:
    temp_directory <- "~/temp" #tempdir()
    dir.create(temp_directory)
    for (dir in folders_zip_data){
      x <- tail(unlist(strsplit(dir, "/")), n=1)
      folder_name <- unlist(strsplit(x, "[.]zip"))
      x <- unzip(dir, exdir=paste(temp_directory, folder_name, sep="/"))
    }
    files_bruker <- list.files(temp_directory, include.dirs=T, full.names=T, recursive=T)
  }
  else files_bruker <- files_study 
  directories_to_read <- grep("/pdata/1$", files_bruker, value=T)
  
  #Get samples names:
  samples_dirs_list <- c() #Name of the directories of each sample
  samples_names_list <- c() #Name of each sample
  #Sample name in index 1 of samples_names_list corresponds to sample directory in index 1 of samples_dirs_list ^
  
  if (!is.null(samples.names)){
    samples_names_df <- read.csv(samples.names, header=F, stringsAsFactors=F)
    if(length(grep("zip$", samples_names_df[[2]], value=T))!=0){
      samples_names_df[[2]] <- substr(samples_names_df[[2]], 1, nchar(as.character(samples_names_df[[2]]))-4)
    }
    samples_dirs_list <- samples_names_df[,2]
    samples_names_list <- samples_names_df[,1]
  }
  else{
    for(dir in directories_to_read){
      n <- tail(strsplit(dir, "/")[[1]], n=4)
      n <- n[1]
      samples_names_list <- c(samples_names_list, n)
      samples_dirs_list <- c(samples_dirs_list, n)
    }
  }
  if (length(directories_to_read) != length(samples_names_list)) stop("File structure is incorrect. Maybe you have more than one spectrum in a sample folder or the samples folders inside each .zip have repeated names. Please look at the documentation.")
  
  #Read metadata file, if given:
  if(!is.null(metadata_file)){
    cat("Reading Metadata file\n")
    metadata <- specmine::read_metadata(metadata_file, header.col=m.header_col, header.row=m.header_row, sep=m.sep)
  }
  else metadata <- NULL
  
  
  spectra.list <- list()
  samples <- c()
  #Read each spectra files:
  for (i in 1:length(directories_to_read)){
    dir <- directories_to_read[i]
    sample_name <- samples_names_list[match(tail(unlist(strsplit(dir, "/")), n=4)[1], samples_dirs_list)]
    
    samples <- c(samples,sample_name)
    
    cat("Reading sample ", sample_name, " in ", dir)
    cat("\n")
    
    sample <- mrbin::readBruker(folder=dir,dimension="2D")
    matrix_sample <- sample$currentSpectrum
    if (!is.null(matrix_sample)){
      rownames(matrix_sample) <- round(as.numeric(rownames(matrix_sample)),digits=2)
      colnames(matrix_sample) <- round(as.numeric(colnames(matrix_sample)),digits=2)
      
      
      spectra.list[[sample_name]] <- matrix_sample
    }
  }
  
  cat("Creating dataset (this may take a while)\n")
  
  dataset <- create_2d_dataset(spectra.list, type="2d-nmr", metadata = metadata, description = description, label.x = label.x,
                               label.y = label.y, label.values = label.values, F1 = as.numeric(rownames(spectra.list[[1]])), F2 = as.numeric(colnames(spectra.list[[1]])))
  
  
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

read_varian_spectrum_raw=function(directory, zero_filling=TRUE, apodization=TRUE){
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package reticulate needed for this function to work. Please install it: install.packages('reticulate')",
         call. = FALSE)
  }
  
  files=list.files(directory, recursive=TRUE, full.names = TRUE)
  fid_file=grep(".*fid$", files, value=TRUE)
  procpar_file=grep(".*procpar.*$", files, value=TRUE)
  
  if(!reticulate::py_available()){
    p=reticulate::py_discover_config()
    python_versions=p$python_versions
    python_path=""
    if(is_windows() & length(python_versions)>0) python_path=python_versions[length(python_versions)]
    else python_path=grep(".*python3.*", p$python_versions, value=TRUE)
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
                                 metadata_file=NULL, m.header_col=TRUE, m.header_row=TRUE, m.sep=",",
                                 samples.names=NULL, zero_filling=TRUE, apodization=TRUE, zipped=TRUE,
                                 description="", label.x="ppm", label.values="intensity"){
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package reticulate needed for this function to work. Please install it: install.packages('reticulate')",
         call. = FALSE)
  }
  
  spectra.list=list()
  
  #Get data paths of the files to read:
  if (zipped){
    files_study=list.files(varian_spectra_directory, include.dirs=TRUE, full.names=TRUE)
    folders_zip_data=grep("fid.zip$", files_study, value=TRUE)
    #Read each spectra files:
    temp_directory="~/temp" #tempdir()
    dir.create(temp_directory)
    for (dir in folders_zip_data){
      unzip(dir, exdir=temp_directory)
    }
    directories_to_read=grep(".*[.]fid", list.dirs(temp_directory, full.names=TRUE), value=TRUE)
  }
  else directories_to_read=grep(".*[.]fid", list.dirs(varian_spectra_directory, full.names=TRUE), value=TRUE)
  
  #Get samples names:
  if (!is.null(samples.names)){
    samples_names=read.csv(samples.names, header=FALSE)
    if(length(grep("zip$", samples_names[[2]], value=TRUE))!=0){
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
  message("Reading files:\n")
  data_matrix_column_names=c()
  for (i in 1:length(directories_to_read)){
    dir=directories_to_read[i]
    message(dir)
    message("\n")
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
  
  if(zipped) unlink(temp_directory, recursive=TRUE)
  
  return(dataset)
  
}

#################################################
#####READ 2D NMR SPECTRA FILES IN VARIAN FORMAT#####
#################################################
##One Varian directory corresponds to one 2D varian spectrum:
###spectrum_directory.fid/
#####-fid
#####-procpar
#####-other folders
#####-other files

read_varian_2dspectrum_raw=function(directory, zero_filling=T, apodization=T){
  
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
  
  
  varian_2draw_py <- system.file("read_varian_2dspec_raw.py", package="specmine")
  env <- new.env()
  reticulate::source_python(varian_2draw_py,  envir = env)
  

  spec <- env$read_varian_spec2d_raw(directory, fid_file, procpar_file, zero_filling=zero_filling, apodization=apodization)
  names(spec) = c("ppms","intensities")
  names(spec$ppms) = c("F2","F1")
  return(spec)
}


read_varian_2dspectra_raw=function(varian_spectra_directory,
                                   metadata_file=NULL, m.header_col=T, m.header_row=T, m.sep=",",
                                   samples.names=NULL, zero_filling=T, apodization=T, zipped=T,
                                   description="", label.x="ppm", label.y="ppm", label.values="intensity"){
  
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
  
  #Read metadata file, if given:
  if(!is.null(metadata_file)){
    cat("Reading Metadata file\n")
    metadata=specmine::read_metadata(metadata_file, header.col=m.header_col, header.row=m.header_row, sep=m.sep)
  }
  else metadata=NULL
  
  #Get spectra, processed accordingly:
  cat("Reading files:\n")
  data_matrix_column_names=c()
  for (i in 1:length(directories_to_read)){
    dir=directories_to_read[i]
    cat(dir)
    cat("\n")
    if (!is.null(samples.names)) sample_name=samples_names[1][samples_names[2]==tail(strsplit(dir, "/")[[1]], n=1)]
    else sample_name=samples_names[i]
    spectrum=read_varian_2dspectrum_raw(dir, zero_filling=zero_filling, apodization=apodization)
    data_2d = as.matrix(spectrum$intensities)
    rownames(data_2d) = round(spectrum$ppms$'F1',digits = 2)
    colnames(data_2d) = round(spectrum$ppms$'F2',digits = 2)
    spectra.list[[sample_name]]= data_2d
    data_matrix_column_names=c(data_matrix_column_names, sample_name)
  }
  
  cat("Creating dataset (this may take a while)\n")
  
  
  dataset=create_2d_dataset(spectra.list, type="2d-nmr", metadata=metadata, description=description, label.x=label.x, label.y=label.y,label.values=label.values, F1 = as.numeric(rownames(spectra.list[[1]])), F2 = as.numeric(colnames(spectra.list[[1]])))
  
  if(zipped) unlink(temp_directory, recursive=T)
  
  return(dataset)
  
}





##########################
#####DETECT NMR PEAKS#####
##########################

detect_nmr_peaks=function(spectrum.list, baseline_treshold=50000){
  
  res=list()
  
  x=t(as.matrix(spectrum.list$intensity))
  peaks_detected_idx=speaq::detectSpecPeaks(x, verbose=FALSE, baselineThresh=baseline_treshold)
  
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
    
    message("Spectrum ", spectrum, " has ", length(detected_peaks$ppm), " peaks.")
    message("\n")
    
    detected_peaks_list[[spectrum]]=as.data.frame(detected_peaks)
    detected_peaks_list[[spectrum]]$ppm=round(detected_peaks_list[[spectrum]]$ppm, 2)
    
  }
  
  final_dataset=specmine::group_peaks(detected_peaks_list, type="nmr-peaks", metadata=dataset$metadata,
                                      method=ap.method, samp.classes=ap.samp.classes, step=ap.step,
                                      label.x=dataset$labels$x, label.values=dataset$labels$val, description=dataset$description)
  
  return(final_dataset)
  
}

#############################
#####DETECT 2D NMR PEAKS#####
#############################

## Internal 2D peak picking function, retrieved from rNMR source code
## Finds points in a matrix that are larger than all surrounding points
## x  -  A numeric matrix containing the range of data to be peak picked
## thresh    - Numeric value specifying the minimum level to be included
## noiseFilt - Integer argument that can be set to 0, 1 or 2; 
##              0 does not apply a noise filter, 1 applies a mild filter
##              (adjacent points in the direct dimension must be above the 
##              noise threshold), 2 applies a strong filter (all adjacent points
##              must be above the noise threshold
## Returns a vector of points defining the local maxima

localMax <- function(x, thresh, noiseFilt) {
  
  nC <- ncol(x)
  nR <- nrow(x)
  if( noiseFilt == 2 )
    x[x < thresh] <- NA  
  
  ## Find row/column local maxes
  if( noiseFilt == 1 ){
    y <- x
    y[ y < thresh ] <- NA
    vMax <- intersect(which(c(NA, y) < c(y, NA)), which(c(NA, y) > c(y, NA))-1)
  }else
    vMax <- intersect(which(c(NA, x) < c(x, NA)), which(c(NA, x) > c(x, NA))-1)
  x <- t(x)
  hMax <- intersect(which(c(NA, x) < c(x, NA)), which(c(NA, x) > c(x, NA))-1)-1
  hMax <- (hMax %% nC * nR) + hMax %/% nC + 1
  
  
  ## Find diagonal maxima
  x <- t(x)
  hvMax <- intersect(vMax, hMax)
  if( noiseFilt == 0 )
    hvMax <- hvMax[ x[hvMax] > thresh ]
  dMax <- cbind(hvMax, hvMax - nR + 1, hvMax - nR - 1, hvMax + nR + 1, 
                hvMax + nR - 1)
  dMax[dMax < 1 | dMax > nC*nR ] <- NA
  dMax <- which(max.col(cbind( x[dMax[,1]], x[dMax[,2]], x[dMax[,3]], 
                               x[dMax[,4]], x[dMax[,5]])) == 1)
  
  return(hvMax[dMax])
}



#Function to change NA values in a data matrix according to a peak list
peaks_to_dataset <- function(empty_data,peaklst,reference) {
  for (i in 1:nrow(peaklst)) {
    if (!(is.na(peaklst$rows[i]) | is.na(peaklst$cols[i]))) {
      empty_data[peaklst$rows[i],peaklst$cols[i]] <- reference[peaklst$rows[i],peaklst$cols[i]]
    }
  }
  empty_data
}



## Function to give a set of peaklists found in 1 spectrum
## Finds pairs of ppms that have a peak
## spectrum  - Numeric matrix. A 2D NMR spectrum.
## threshold - Numeric value. Option to user establish a threshold defining a minimum value to be detected
## noise - Integer argument that can be set to 0, 1 or 2; 
##              0 does not apply a noise filter, 1 applies a mild filter
##              (adjacent points in the direct dimension must be above the 
##              noise threshold), 2 applies a strong filter (all adjacent points
##              must be above the noise threshold
## Returns a data frame containing the pairs (row/column) and the intensity of the peaks detected
peaklist <- function(spectrum, threshold = NULL, noise=0) {
  
  if (is.null(spectrum)) {
    stop("Specrum not provided")
  }
  
  if (!is.null(threshold)) {
    
    #Calculate local Maxes
    result <- localMax(spectrum,thresh = threshold,noiseFilt = noise)
    
    #Translates them to points in dataset
    points <- spectrum[result]
    
    #Associate its point to its row and column
    locations <- as.matrix(lapply(points,function(x)(which(x==spectrum,arr.ind = T))))
    
    #Return the data frame with peak information
    peaks <- data.frame(rows = unlist(lapply(locations,function(x)rownames(spectrum)[x[1,][1]])),cols = unlist(lapply(locations,function(x)(colnames(spectrum)[x[1,][2]]))))
    
  } else {
    #No threshold given so, mean of intensity values used
    
    #SNR calculated
    #threshold <- snr_spectra(spectrum)
    #threshold2 <- max(spectrum)/(mean(spectrum[which(apply(spectrum,1,sd)<mean(spectrum)),]))
    threshold <- mean(base::Filter(isPositive,spectrum))
    
    #Now that threshold is calculated, detect peaks
    #Calculate local Maxes
    result <- localMax(spectrum,thresh = threshold,noiseFilt = noise)
    
    #Translates them to points in dataset
    points <- spectrum[result]
    
    #Associate its point to its row and column
    locations <- as.matrix(sapply(points,function(x)(which(x==spectrum,arr.ind = T))))
    
    sample_rows <- sapply(locations[1,],function(x)rownames(spectrum)[x])
    sample_cols <- sapply(locations[2,],function(x)colnames(spectrum)[x])
    
    
    #Return the data frame with peak information
    peaks <- data.frame(rows = sample_rows, cols = sample_cols)
    
  }
}



## Function to detect peaks in a specmine dataset of 2D NMR data
## Finds peaks across sample, reducing dimensionality
## specmine_2d_dataset  - A list of variables containing at least a 2D matrix for all samples
## thresh - Numeric value. Option to user establish a threshold defining a minimum value to be detected
## noiseFilt - Integer argument that can be set to 0, 1 or 2; 
##              0 does not apply a noise filter, 1 applies a mild filter
##              (adjacent points in the direct dimension must be above the 
##              noise threshold), 2 applies a strong filter (all adjacent points
##              must be above the noise threshold
## negatives - Boolean value to decide if negative ppm values should be considered or not  
## Returns a specmine dataset with only the variables that was found a peak for, a normal 1D

peak_detection2d <- function(specmine_2d_dataset, baseline_thresh=NULL, noiseFilt=0, negatives=F) {
  
  #list of 2D spectra | 1 sample = 1 spectra  
  data <- specmine_2d_dataset$data
  
  #New list to store the new peak detected datasets
  res_data <- list()
  
  for (i in 1:length(data)) {
    sample <- names(data)[[i]]
    
    rows <- nrow(data[[i]])
    cols <- ncol(data[[i]])
    
    #create new matrix filled with NA
    res_data[[sample]] <- matrix(data = rep(NA,rows*cols),nrow=rows,ncol=cols,dimnames = list(rownames(data[[i]]),colnames(data[[i]])))
    
    #Calculate peaks
    peaklist <- peaklist(spectrum = data[[i]],threshold = baseline_thresh, noise = noiseFilt)
    cat(paste('Sample:', sample, 'has', nrow(peaklist),'peaks \n'))
    
    res_data[[sample]] <- peaks_to_dataset(empty_data = res_data[[sample]],peaklst = peaklist,reference = data[[i]])
    
  }
  
  dim_example <- dim(res_data[[1]])
  logical <- unlist(lapply(lapply(res_data,dim),function(x)identical(x,dim_example)))
  if (sum(logical) == length(res_data)) {
    res_data <- simplify2array(res_data)
  } else {
    res_data <- narray::stack(res_data)
  }
  
  
  if (!negatives) {
    col <- grep("-", colnames(res_data))
    row <- grep("-", rownames(res_data))
    #Remove negative columns
    if (length(col)>0) {
      res_data <- res_data[,-(col),]
    }
    #Remove negative rows
    if (length(row)>0) {
      res_data <- res_data[-(row),,]
    }
  }
  
  #Transfrom into a 2D matrix
  dimnames1 <- as.vector(t(outer(dimnames(res_data)[[1]], dimnames(res_data)[[2]], FUN = paste, sep=".")))
  data_2d <- matrix(res_data,prod(dim(res_data)[1:2]),dim(res_data)[3],dimnames = list(dimnames1,dimnames(res_data)[[3]]))
  # data_2d <- t(two.d.array(res_data))
  
  #Make unique variable names due to repeated ppm values
  rownames(data_2d) <- base::make.names(rownames(data_2d),unique = TRUE)
  
  #Find which rows have a full NA vector, in order to remove it
  indexes <- which(rowSums(is.na(data_2d)) == ncol(data_2d))
  
  
  #Filter those rows out
  data_2d <- data_2d[-indexes,]
  
  
  #Create a new list where the data has only rows containing peaks
  dataset <- specmine::create_dataset(data_2d, type ='nmr-peaks', metadata = specmine_2d_dataset$metadata, description = specmine_2d_dataset$description,
                                      label.x = 'F1 x F2 ppm', label.values = 'intensity', sample.names = names(specmine_2d_dataset$data))
  
  
  
  return(dataset)
  
}


