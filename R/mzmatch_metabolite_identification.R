
mzmatch.identify.metabolites = function(ionisation="detect", data.folder = NULL, xSet = NULL, backend = "Ramp",
										adducts = "M+H,M+ACN+Na,M+Na,M+K,M+ACN+H"){
										
	require(mzmatch.R)
	mzmatch.init(version.1=FALSE)
	
	if (is.null(xSet)){
		require(xcms)
		files = list.files(data.folder, recursive=T, full.names=TRUE)
		xSet <- xcmsSet(files, method='centWave', ppm=2, peakwidth=c(5,100), 
                  snthresh=5, prefilter=c(3,1000), integrate=1, mzdiff=0.001, 
                  verbose.columns=TRUE, fitgauss=FALSE, nSlaves=4)
	}
	
	peakml.files = gsub("(.+)/([^/]+\\..+$)","\\2", xSet@filepaths)
	peakml.files = paste("./peakml/", peakml.files, sep="")
	
	PeakML.xcms.write.SingleMeasurement (xset=xSet,filename=peakml.files,
                                       ionisation=ionisation,addscans=2,
                                       writeRejected=FALSE,ApodisationFilter=TRUE)
    
	mzmatch.ipeak.Combine(i=paste(peakml.files,collapse=","),v=T,rtwindow=30,
						  o="./mzmatch/combined.peakml",combination="set",ppm=5)
	mzmatch.ipeak.filter.NoiseFilter (i="./mzmatch/combined.peakml",o="./mzmatch/combined_noisef.peakml",
									  v=T,codadw=0.8)
	mzmatch.ipeak.filter.SimpleFilter(i="./mzmatch/combined_noisef.peakml", 
									  o="./mzmatch/combined_sfdet.peakml", mindetections=3)
	mzmatch.ipeak.filter.SimpleFilter(i="./mzmatch/combined_sfdet.peakml", 
									  o="./mzmatch/combined_highintensity.peakml", 
									  minintensity=100000)
	PeakML.GapFiller(filename = "./mzmatch/combined_highintensity.peakml", ionisation = "detect", 
					 outputfile = "./mzmatch/highintensity_gapfilled.peakml", 
					 ppm=5, rtwin = 0, backend = backend)
	mzmatch.ipeak.sort.RelatedPeaks (i="./mzmatch/highintensity_gapfilled.peakml",v=T,
									 o="./mzmatch/mzMatch_output.peakml",
									 basepeaks="./mzmatch/mzMatch_basepeaks.peakml",ppm=3,
									 rtwindow=6)
									 
	annot <- paste("relation.id,relation.ship,codadw,charge")
	mzmatch.ipeak.convert.ConvertToText (i="./mzmatch/mzMatch_output.peakml",
										 o="./mzmatch/mzMATCHoutput.txt",v=T,annotations=annot)

	DBS <- dir(paste(find.package("mzmatch.R"), "/dbs", sep=""),
			   full.names=TRUE)
	DBS <- paste(DBS,collapse='","')
	DBS <- paste('"',DBS,'"',sep="")

	mzmatch.ipeak.util.Identify(i="./mzmatch/mzMatch_output.peakml", v=T,
								o="./mzmatch/final_combined_related_identified.peakml", ppm=3, 
								databases=DBS, adducts = adducts)


	mzmatch.ipeak.convert.ConvertToText (
	  i="./mzmatch/final_combined_related_identified.peakml",
	  o= "./mzmatch/final_combined_related_identified.txt", databases=DBS,
	  annotations="identification,ppm,adduct,relation.ship")
    

}
