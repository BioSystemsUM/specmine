# A commonly used value is 30 (seconds) for LC-MS and 4 (seconds) for GC-MS spectra (fwhm).
# read LC/GC-MS spectra(.netCDF, .mzXML, mzData)
# use functions in XCMS package
read.ms.spec<-function(folder.name, metadata = NULL, profmethod='bin', fwhm=30, bw=30){
  suppressMessages(require(xcms));
  files <- list.files(folder.name, recursive=T, full.names=TRUE);
  if (is.null(metadata)) xset <- xcmsSet(files, profmethod = "bin", fwhm=fwhm);
  else xset = xcmsSet(files, profmethod = "bin", fwhm = fwhm, sclass = metadata[,1);
  xset<-group(xset, bw=bw);
  xset
}



# retention time correction for LC/GC-MS spectra
ms.rt.correction<-function(xset, bw=30){
  xset2<-retcor(xset)
  # re-group peaks after retention time correction
  xset2<-group(xset2, bw=30)
  xset2
}



# fill in missing peaks
ms.fill.peaks<-function(xset){
  xset2<-fillPeaks(xset);
  xset2
}


# into:  integrated area of original (raw) peak
# intf:  integrated area of filtered peak
# maxo:  maximum intensity of original (raw) peak
# maxf:  maximum intensity of filtered peak

ms.create.matrix<-function(xset, intvalue = "into"){
  values <- groupval(xset, "medret", value = intvalue);
  values
}

# A commonly used value is 30 (seconds) for LC-MS and 4 (seconds) for GC-MS spectra (fwhm).
read.ms.spectra = function(folder.name, type = "undefined", filename.meta= NULL, description = "", prof.method='bin', fwhm=30, bw=30, intvalue = "into", header.col.meta = TRUE, header.row.meta = TRUE, sep.meta = ","){
	if (!is.null(filename.meta))
		metadata = read.metadata(filename.meta, header.col = header.col.meta, header.row = header.row.meta, sep = sep.meta)
	else metadata = NULL
	xset = read.ms.spec(folder.name, metadata = metadata, prof.method, fwhm, bw)
	xset = ms.rt.correction(xset, bw)
	xset = ms.fill.peaks(xset)
	mat = ms.create.matrix(xset, intvalue)
	dataset = create.dataset(mat, type = type, metadata = metadata, description = description, 
							 label.x = "mz/rt", label.values = "intensity", xSet = xset)
	dataset
  
}


