# A commonly used value is 30 (seconds) for LC-MS and 4 (seconds) for GC-MS spectra (fwhm).
# read LC/GC-MS spectra(.netCDF, .mzXML, mzData)
# use functions in XCMS package
read.ms.spec<-function(folder.name, profmethod='bin', fwhm=30, bw=30){
  suppressMessages(require(xcms));
  files <- list.files(folder, recursive=T, full.names=TRUE);
    
  xset <- xcmsSet(msfiles, profmethod = "bin", fwhm=fwhm);
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
read.ms.spectra = function(folder.name, type = "undefined", metadata = NULL, description = "", prof.method='bin', fwhm=30, bw=30, intvalue = "into"){
	xset = read.ms.spec(folder.name, prof.method, fwhm, bw)
	xset = ms.rt.correction(xset, bw)
	xset = ms.fill.peaks(xset)
	mat = ms.create.matrix(xset, intvalue)
	dataset = create.dataset(mat, type = type, metadata = metadata, description = description, 
							 label.x = "mz/rt", label.values = "intensity")
	dataset
  
}


