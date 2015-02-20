"read.dataset.spc" = function(folder.data, filename.meta= NULL, type = "undefined", description = "", modified = F,
                              label.x = NULL, label.values = NULL,  
                              header.col.meta = TRUE, header.row.meta = TRUE, sep.meta = ","){
								
	if (!is.null(filename.meta))
		metadata = read.metadata(filename.meta, header.col = header.col.meta, header.row = header.row.meta, sep = sep.meta)
	else metadata = NULL
	
	data.spc = read.data.spc(folder.data, modified = modified)
	
	freqs = data.spc[[1]]$wavelength # get frequencies from first spectrum
	datamat = matrix(data = NA, nrow = length(freqs), ncol = length(data.spc))

	for (i in 1:length(data.spc)) datamat[,i] = data.spc[[i]]$spc

	rownames(datamat) = as.character(freqs)
	colnames(datamat) = names(data.spc)
	
	if (is.null(label.x)) label.x = data.spc[[1]]$labels$.wavelength
	if (is.null(label.values)) label.values = data.spc[[1]]$labels$spc

	dataset = create.dataset(datamat, type = type, metadata = metadata, description = description, 
						   label.x = label.x, label.values = label.values)
	dataset
}

"read.data.spc" = function(foldername, modified = F)
{
  require(hyperSpec)
  filenames = dir(foldername, pattern=".[Ss][Pp][Cc]$", full.name=TRUE)
  sampleList = list()
  sampleNames = c()
  snames <- gsub("\\.[^.]*$", "", basename(filenames));
  for (i in 1:length(filenames)) {
    print(paste("Reading sample ", filenames[i]))
    if (!modified){
		sampleList[[i]] = read.spc(filenames[i], no.object = T)
	} else {
		sampleList[[i]] = read.spc.modified(filenames[i], no.object = T)
	}
  }
  sampleNames = snames
  names(sampleList) = sampleNames
  sampleList
}

get.samples.names.spc = function(foldername){
  files = list.files(foldername,pattern=".[Ss][Pp][Cc]$", recursive = TRUE, full.name= TRUE)
  samples.names = gsub("\\.[^.]*$", "",basename(files))
  samples.names
}
