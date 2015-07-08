
mzmatch.identify.metabolites = function(ionisation="detect", data.folder = NULL, xSet = NULL, backend = "Ramp",
										adducts = "M+H,M+ACN+Na,M+Na,M+K,M+ACN+H"){
										
	require(mzmatch.R)
	require(rJava)
	mzmatch.init(version.1=FALSE)
	
	if (is.null(xSet)){
		require(xcms)
		files = list.files(data.folder, recursive=T, full.names=TRUE)
		xSet <- xcmsSet(files, method='centWave', ppm=2, peakwidth=c(5,100), 
                  snthresh=5, prefilter=c(3,1000), integrate=1, mzdiff=0.001, 
                  verbose.columns=TRUE, fitgauss=FALSE, nSlaves=4)
	}
	
	peakml.files = gsub("(.+)/([^/]+\\..+$)","\\2", xSet@filepaths)
	peakml.files = gsub("\\..+$", "\\.peakml", peakml.files)
	peakml.files = paste("./peakml/", peakml.files, sep="")
	dir.create("peakml")
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
	PeakML.GapFiller.modified(filename = "./mzmatch/combined_highintensity.peakml", ionisation = ionisation, 
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
    
	metabolites = PeakML.Read("./mzmatch/final_combined_related_identified.peakml")
	metabolites
}


PeakML.GapFiller.modified = function (filename, ionisation = "detect", Rawpath = NULL, outputfile, 
    ppm = 0, rtwin = 0, nSlaves = 1, fillAll = FALSE, backend = "Ramp") 
{
    version.1 <- get("version.1", envir = .GlobalEnv)
    FillinPeaks <- function(peaknum) {
        whichpeakset <- numchromsexpected[fillinnums[peaknum], 
            1]
        subtable <- PeakMLdata$peakDataMtx[PeakMLdata$peakDataMtx[, 
            10] == whichpeakset, ]
        subtable <- rbind(subtable, NULL)
        rt_start <- min(subtable[, 5]) - rtwin
        rt_finis <- max(subtable[, 6]) + rtwin
        if (rt_finis > max(correctedRT)) {
            rt_finis <- max(correctedRT)
        }
        if (rt_start > max(correctedRT)) {
            rt_start <- max(correctedRT)
        }
        if (rt_finis < min(correctedRT)) {
            rt_finis <- min(correctedRT)
        }
        if (rt_start < min(correctedRT)) {
            rt_start <- min(correctedRT)
        }
        mz_start <- min(subtable[, 2])
        mz_finis <- max(subtable[, 3])
        mz_start <- mz_start - (mz_start * ppm/10^6)
        mz_finis <- mz_finis + (mz_finis * ppm/10^6)
        scan_start <- which(correctedRT >= rt_start)[1] - 1
        if (scan_start == 0) 
            scan_start = 1
        scan_finis <- which(correctedRT >= rt_finis)[1]
        C <- try(PeakML.Methods.getRawMat(allRawPeaks, scan_start, 
            scan_finis, mz_start, mz_finis, correctedRT, uncorrectedRT))
        if (class(C) == "try-error") {
            C <- c(1, 1, 1, 1, 1)
        }
        C <- rbind(C, NULL)
        if (nrow(C) <= 3 | length(unique(C[, 5])) <= 3) {
            scanids <- c(-1, -1, -1)
            retentiontimes <- c(-1, -1, -1)
            masses <- c(-1, -1, -1)
            intensities <- c(-1, -1, -1)
        }
        else {
            scanids <- C[, 3]
            retentiontimes <- C[, 2]
            masses <- C[, 4]
            intensities <- C[, 5]
        }
        OUT <- rbind(masses, intensities, retentiontimes, scanids - 
            1)
    }
    st <- system.time(PeakMLdata <- PeakML.Read(filename, ionisation, 
        Rawpath))
    ionisation <- PeakMLdata$massCorrection[[2]]
    massCorrection <- PeakMLdata$massCorrection[[1]]
    samplenames <- PeakMLdata$sampleNames
    rawdatafullpaths <- PeakMLdata$rawDataFullPaths
    if (is.null(rawdatafullpaths)) {
        cat("Some of the raw data files are not accessible, we will not be able to fill in missing peaks. Please set \"Rawpath\" argument with location where files can be located\n")
        stop()
    }
    numchromsexpected <- unlist(lapply(1:max(PeakMLdata$peakDataMtx[, 
        10]), function(x) rep(x, length(samplenames))))
    numchromsexpected <- cbind(numchromsexpected, NA, NA, NA)
    for (setnum in 1:max(PeakMLdata$peakDataMtx[, 10])) {
        inset <- c(1:length(samplenames))
        rownums <- which(PeakMLdata$peakDataMtx[, 10] == setnum)
        hit <- PeakMLdata$peakDataMtx[rownums, 9]
        oneEqualsSetnum = which(numchromsexpected[, 1] == setnum)
        numchromsexpected[oneEqualsSetnum, 2] <- as.numeric(inset %in% 
            hit)
        missed <- which(inset %in% hit == FALSE)
        if (length(missed) > 0) {
            detectedpeaks <- c(rep(1, length(hit)), rep(0, length(missed)))
            hit <- append(hit, missed)
            rownums <- append(rownums, rep(0, length(missed)))
        }
        else {
            detectedpeaks <- rep(1, length(hit))
        }
        numchromsexpected[oneEqualsSetnum, 3] <- hit
        numchromsexpected[oneEqualsSetnum, 2] <- detectedpeaks
        numchromsexpected[oneEqualsSetnum, 4] <- rownums
    }
    colnames(numchromsexpected) <- NULL
    chromslist <- vector("list", nrow(numchromsexpected))
    if (fillAll == TRUE) {
        detectedchromatograms <- numchromsexpected[, 2]
        numchromsexpected[, 2] <- 0
    }
    notdetected <- which(numchromsexpected[, 2] == 0)
    whichfiles <- numchromsexpected[notdetected, 3]
    samplenums <- unique(whichfiles)
    if (length(samplenums != 0)) {
        for (filenum in 1:length(samplenums)) {
            samplefile <- samplenums[filenum]
            cat("Working on file: ", rawdatafullpaths[samplefile], 
                "\n")
            rawfile <- openMSfile(rawdatafullpaths[samplefile], 
                verbose = FALSE, backend = backend)
            allRawPeaks <- peaks(rawfile)
            correctedRT <- as.numeric(PeakMLdata$correctedRTList[[samplefile]])
            uncorrectedRT <- header(rawfile)$retentionTime
            if (all(correctedRT == uncorrectedRT)) {
                rtCorrection <- FALSE
            }
            else {
                rtCorrection <- TRUE
            }
            fillinnums <- notdetected[whichfiles == samplefile]
            isSnow <- FALSE
            if (nSlaves > 1) {
                HIT <- grep("snow", installed.packages()[, 1])
                if (!is.null(HIT)) {
                  isSnow == TRUE
                }
                else {
                  cat("Pleae install package snow to use multiple processors. \n We will continue with a single processor for the time being.", 
                    "\\n")
                }
            }
            if (isSnow == TRUE) {
                if (filenum == 1) {
                  cat("Package snow loaded.", "\n")
                }
                cl <- makeCluster(nSlaves)
                assign("rtwin", rtwin, envir = .GlobalEnv)
                assign("rawfile", rawfile, envir = .GlobalEnv)
                assign("numchromsexpected", numchromsexpected, 
                  envir = .GlobalEnv)
                assign("fillinnums", fillinnums, envir = .GlobalEnv)
                assign("PeakMLdata$peakDataMtx", PeakMLdata$peakDataMtx, 
                  envir = .GlobalEnv)
                assign("ppm", ppm, envir = .GlobalEnv)
                assign("FillinPeaks", FillinPeaks, ppm, envir = .GlobalEnv)
                assign("rawMat", rawMat, envir = .GlobalEnv)
                clusterExport(cl, list = c("rtwin", "rawfile", 
                  "numchromsexpected", "fillinnums", "PeakMLdata$peakDataMtx", 
                  "ppm", "FillinPeaks", "rawMat"))
                filledlist <- parLapply(cl, c(1:length(fillinnums)), 
                  FillinPeaks)
                stopCluster(cl)
            }
            else {
                filledlist <- lapply(1:length(fillinnums), FillinPeaks)
            }
            for (i in 1:length(filledlist)) {
                chromslist[[fillinnums[i]]] <- filledlist[[i]]
            }
            mzR::close(rawfile)
            rm(filledlist)
        }
    }
    project <- .jnew("peakml/util/rjava/Project", samplenames, 
        rawdatafullpaths, as.character(PeakMLdata$phenoData))
    .jcall(project, returnSig = "V", method = "addHeaderAnnotation", 
        as.character("peakproc"), as.character("XCMS_Gapfilled"))
    for (measurementid in 1:length(samplenames)) {
        for (scannum in 1:length(PeakMLdata$correctedRTList[[measurementid]])) {
            .jcall(project, returnSig = "V", method = "addScanInfo", 
                as.integer(measurementid - 1), as.numeric(PeakMLdata$correctedRTList[[measurementid]][scannum]), 
                as.character(ionisation))
            .jcall(project, returnSig = "V", method = "addScanAnnotation", 
                as.integer(measurementid - 1), as.integer(scannum - 
                  1), as.character("RT_raw"), as.character(PeakMLdata$rawRTList[[measurementid]][scannum]))
        }
    }
    for (i in 1:length(chromslist)) {
        if (numchromsexpected[i, 2] == 0) {
            chrom <- chromslist[[i]]
        }
        else {
            ind <- numchromsexpected[i, 4]
            chrom <- PeakMLdata$chromDataList[[ind]]
        }
        .jcall(project, returnSig = "V", method = "addMassChromatogram", 
            as.integer(numchromsexpected[i, 3] - 1), as.integer(chrom[4, 
                ]), as.numeric(chrom[3, ]), as.numeric(chrom[1, 
                ]), as.numeric(chrom[2, ]), as.character(ionisation))
    }
    setindexes <- vector("list", length(unique(numchromsexpected[, 
        1])))
    for (indexnumber in 1:length(setindexes)) {
        setindexes[[indexnumber]] <- which(numchromsexpected[, 
            1] == indexnumber)
    }
    for (ind in 1:length(setindexes)) {
        .jcall(project, returnSig = "V", method = "addPeakSet", 
            as.integer(setindexes[[ind]] - 1))
    }
    if (!is.null(PeakMLdata$GroupAnnotations)) {
        PeakML.Methods.writeGroupAnnotations(project, PeakMLdata$GroupAnnotations)
    }
    .jcall(project, returnSig = "V", method = "write", outputfile)
}
