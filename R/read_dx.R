# functions to read JDX spectra files 

"read_dataset_dx" = function(folder.data, filename.meta= NULL, type = "undefined", description = "", 
                              label.x = NULL, label.values = NULL,  
                              header.col.meta = TRUE, header.row.meta = TRUE, sep.meta = ",")
{
  if (!is.null(filename.meta))
    metadata = read_metadata(filename.meta, header.col = header.col.meta, header.row = header.row.meta, sep = sep.meta)
  else metadata = NULL
  
  data.dx = read_data_dx(folder.data)
  
  freqs = data.dx[[1]][[4]][,1] # get frequencies from first spectrum
  datamat = matrix(data = NA, nrow = length(freqs), ncol = length(data.dx))

  for (i in 1:length(data.dx)) datamat[,i] = data.dx[[i]][[4]][,2]
  
  rownames(datamat) = as.character(freqs)
  colnames(datamat) = names(data.dx)
  
  dataset = create_dataset(datamat, type = type, metadata = metadata, description = description, 
                           label.x = label.x, label.values = label.values)
  dataset
}

"read_data_dx" = function(foldername, debug = 0)
{
  filenames = dir(foldername, pattern=".[Dd][Xx]$", full.names=TRUE)
  sampleList = list()
  sampleNames = c()
  snames <- gsub("\\.[^.]*$", "", basename(filenames));
  for (i in 1:length(filenames)) {
    print(paste("Reading sample ", filenames[i]))
    sampleList[[i]] = readJDX::readJDX(filenames[i], debug = debug, SOFC = FALSE) #If TRUE it fails for files that don't have 'firstY' metadata
  }
  sampleNames = snames
  names(sampleList) = sampleNames
  sampleList
}

"get_samples_names_dx" = function(foldername){
  files = list.files(foldername,pattern=".[Dd][Xx]$", recursive = TRUE, full.names= TRUE)
  samples.names = gsub("\\.[^.]*$", "",basename(files))
  samples.names
}

# "readJDX_specmine" = function (file = "", debug = FALSE) 
# {
#   if (!requireNamespace("gsubfn", quietly = TRUE)) {
#     stop("You need to install package gsubfn to use this option")
#   }
#   if (file == "") stop("No file specified")
#   jdx <- readLines(file)
#   cmpd <- grep("^##TITLE=.*$", jdx)
#   if (debug) cat("CMPD number:" , cmpd)
#   if (cmpd > 1) stop("Compound data sets not supported")
#   if (cmpd == 0) stop("This may not be a JCAMP-DX file")
#   ntup <- grepl("^##NTUPLES", jdx)
#   if (any(ntup)) 
#     stop("This looks like NMR data with real & imaginary parts, which is not supported")
#   if (debug) cat("\nFile = ", file, "\n")
#   spcstart <- grep("^##XYDATA=\\s*\\(X\\+\\+\\(Y\\.\\.Y\\)\\)$",  jdx) 
#   if (debug) cat ("SPCSTART" , spcstart, "\n")
#   if (spcstart == 0) stop("Couldn't find the data block start (see ?readJDX for supported formats)")
#   spcstart = spcstart + 1
#   
#   if (debug) print(jdx[1:(spcstart - 2)])
#   
#   spcend <- grep("^##END=", jdx) - 1
#   if (debug) cat ("SPCEND" , spcend, "\n")
#   if (spcend == 0) stop("Couldn't find the data block end")
# 
#   if (!length(spcstart) == 1L | !length(spcend) == 1L) 
#     stop("Problem with delimiting data block")
#   if (!spcstart < spcend) 
#     stop("End of data block in the wrong place")
# 
#   yValues <- jdx[spcstart:spcend]
#   for (n in 1:length(yValues)) {
#     yValues[n] <- sub("\\s*(\\+|-)*[[:digit:]]+(\\.|,)?[[:digit:]]*\\s*", "", yValues[n])
#   }
#   yValues <- paste(yValues, collapse = " ")
#   yValues <- gsub("\\+", " ", yValues)
#   yValues <- gsub("-", " -", yValues)
#   yValues <- sub("\\s*", "", yValues)
#   yValues <- gsub(",", ".", yValues)
#   yValues <- strsplit(yValues, split = "\\s+")
#   yValues <- as.numeric(unlist(yValues))
#   firstX <- grep("^##FIRSTX=", jdx)
#   if (firstX == 0)stop("Couldn't find FIRSTX")
#   firstX <- jdx[firstX]
#   firstX <- gsubfn::gsubfn("##FIRSTX=", replacement = "", firstX)
#   firstX <- sub(",", ".", firstX)
#   firstX <- as.numeric(firstX)
#   lastX <- grep("^##LASTX=", jdx)
#   if (lastX == 0) stop("Couldn't find LASTX")
#   lastX <- jdx[lastX]
#   lastX <- gsubfn::gsubfn("##LASTX=", replacement = "", lastX)
#   lastX <- sub(",", ".", lastX)
#   lastX <- as.numeric(lastX)
#   npoints <- grep("^##NPOINTS=", jdx)
#   if (npoints == 0) stop("Couldn't find NPOINTS")
#   npoints <- jdx[npoints]
#   npoints <- gsubfn::gsubfn("##NPOINTS=", replacement = "", npoints)
#   npoints <- as.integer(npoints)
#   if (debug) cat("NPOINTS = ", npoints, "\n")
#   if (debug)  cat("Actual no. data points found  = ", length(yValues), "\n")
#   if (!npoints == length(yValues)) stop("NPOINTS and length of data block don't match")
#   if (debug) cat("firstX = ", firstX, "\n")
#   if (debug) cat("lastX = ", lastX, "\n")
#   yFac <- grep("^##YFACTOR=", jdx)
#   if (yFac == 0) stop("Couldn't find YFACTOR")
#   yFac <- gsubfn::gsubfn("##YFACTOR=", replacement = "", jdx[yFac])
#   yFac <- sub(",", ".", yFac)
#   yFac <- as.numeric(yFac)
#   if (debug) cat("yFac = ", yFac, "\n")
#   yValues <- yValues * yFac
#   actDX <- (lastX - firstX)/(npoints - 1)
#   xValues <- seq(firstX, lastX, by = actDX)
#   res <- data.frame(x = xValues, y = yValues)
# }
