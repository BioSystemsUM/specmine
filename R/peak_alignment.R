
# group peaks with one of two methods
# method: "own" - own algorithm, "metaboanalyst" - metaboanalyst algorithm
# samp.classes: column of metadata dataframe - needed in case of "metaboanalyst"
# step used in "own" algorithm
# returns dataset using standard structure
"group.peaks" = function(sample.list, type, method = "own", metadata = NULL, samp.classes= 1, 
                         description = "", label.x = NULL, label.values = NULL, step = 0.03) {
  if (type == "nmr-peaks"){
	mzwid = 0.03
	bw = 10
  } else if (type == "lcms-peaks"){
	mzwid = 0.25
	bw = 30
  } else if (type == "gcms-peaks"){
	mzwid = 0.25
	bw = 5
  } else stop(paste("The ", type," type is not supported!", sep = ""))
	
  if (method == "own"){
		merged.peaks = merge.eq.peaks.samplelist(sample.list)
		samples.df = get.all.intensities(merged.peaks)
		data = group.peaks.own(samples.df, step)
	} 
  else if (method == "metaboanalyst"){
		data = group.peaks.metaboanalyst(sample.list, metadata[,samp.classes], names(sample.list), mzwid = mzwid, bw = bw)
	}
	create.dataset(as.matrix(data), metadata = metadata, type = type, 
                 description = description, label.x = label.x, label.values = label.values)
}


# grouping peaks: our algorithm
"group.peaks.own" = function(samples.df, step = 0.03)
{
  freqs = as.numeric(rownames(samples.df))
  st_intervals = create.intervals(freqs, step)
  initialized = F
  for (cent in st_intervals)
  {
    rows = samples.df[freqs >= cent & freqs < cent + step,]
    if (nrow(rows)==1)
    {
      if (initialized) {
        res = rbind(res, rows)
      }
      else { 
        res = rows
        initialized = T
      }
    }
    else if (nrow(rows) > 1)
    {
      newpeak = c()
      for(i in 1:ncol(samples.df)) {
        if (sum(!is.na(rows[,i])) > 0) newpeak[i] = sum(rows[,i], na.rm=T)
        else newpeak[i] = NA
      }
      if (initialized) {
        res = rbind(res, newpeak)
      }
      else {
        res = data.frame(t(newpeak))
        colnames(res) = colnames(samples.df)
        initialized = T
      }
      fs = c()
      for(k in 1:nrow(rows))
      {
        f = as.numeric(rownames(rows)[k])
        ocs = sum(!is.na(rows[k,]))
        fs = c(fs, rep(f, each = ocs))
      }
      m = round(median(fs),2)
      rownames(res)[nrow(res)] = as.character(m)
    }
  }
  names(res) = names(samples.df)
  res
}

# auxiliary function for our own grouping algorithm

"create.intervals" = function(freqs, step = 0.03)
{
  st = freqs[1]
  res = c(st)
  index = 1
  while(index < length(freqs))
  {
    while(index <= length(freqs) & freqs[index] < st + step) index = index + 1
    if (index <= length(freqs)) {
      st = freqs[index]
      res = c(res, st)
    }
  }
  res
}



# Group peak list based on position using xcms algorithm (align peaks wrt rt and mz)
# NMR peaks change ppm -> mz and add dummy rt
# 2-col MS need to add dummy rt
# 3-col MS can be used directly
# default mzwid MS 0.25 m/z, NMR 0.03 ppm
# bw 30 for LCMS, 5 for GCMS
group.peaks.metaboanalyst<-function(peaklist, samp.classes, samp.names,  mzwid = 0.03, bw = 10, minfrac = 0.5, minsamp = 1, max = 50) {
	require('xcms')
    samples <- samp.names;
    classlabel <- samp.classes;
    classnames <- levels(classlabel)

    classlabel <- as.vector(unclass(classlabel))
    classnum <- integer(max(classlabel))
    for (i in seq(along = classnum)){
        classnum[i] <- sum(classlabel == i)
    }

    peakmatrix = create.metaboanalyst.mat(peaklist)
    porder <- order(peakmatrix[,"ppm"])
    peakmat <- peakmatrix[porder,,drop=F]
    rownames(peakmat) <- NULL
    retrange <- range(peakmat[,"rt"])

    minpeakmat <- min(classnum)/2

    mass <- seq(peakmat[1,"ppm"], peakmat[nrow(peakmat),"ppm"] + mzwid, by = mzwid/2)
    masspos <- findEqualGreaterM(peakmat[,"ppm"], mass)

    groupmat <- matrix(nrow = 512, ncol = 7 + length(classnum))
    groupindex <- vector("list", 512)

    endidx <- 0
    num <- 0
    gcount <- integer(length(classnum))
    for (i in seq(length = length(mass)-2)) {
        startidx <- masspos[i]
        endidx <- masspos[i+2]-1
        if (endidx - startidx + 1 < minpeakmat)
            next
        speakmat <- peakmat[startidx:endidx,,drop=FALSE]
        den <- density(speakmat[,"rt"], bw, from = retrange[1]-3*bw, to = retrange[2]+3*bw)
        maxden <- max(den$y)
        deny <- den$y
        gmat <- matrix(nrow = 5, ncol = 2+length(classnum))
        snum <- 0
        while (deny[maxy <- which.max(deny)] > maxden/20 && snum < max) {
            grange <- xcms:::descendMin(deny, maxy)
            deny[grange[1]:grange[2]] <- 0
            gidx <- which(speakmat[,"rt"] >= den$x[grange[1]] & speakmat[,"rt"] <= den$x[grange[2]])
            gnum <- classlabel[unique(speakmat[gidx,"sample"])]
            for (j in seq(along = gcount))
                gcount[j] <- sum(gnum == j)
            if (! any(gcount >= classnum*minfrac & gcount >= minsamp))
                next
            snum <- snum + 1
            num <- num + 1
            ### Double the size of the output containers if they're full
            if (num > nrow(groupmat)) {
                groupmat <- rbind(groupmat, matrix(nrow = nrow(groupmat), ncol = ncol(groupmat)))
                groupindex <- c(groupindex, vector("list", length(groupindex)))
            }
            groupmat[num, 1] <- median(speakmat[gidx, "ppm"])
            groupmat[num, 2:3] <- range(speakmat[gidx, "ppm"])
            groupmat[num, 4] <- median(speakmat[gidx, "rt"])
            groupmat[num, 5:6] <- range(speakmat[gidx, "rt"])
            groupmat[num, 7] <- length(gidx)
            groupmat[num, 7+seq(along = gcount)] <- gcount
            groupindex[[num]] <- sort(porder[(startidx:endidx)[gidx]])
        }
    }
    colnames(groupmat) <- c("ppmmed", "ppmmin", "ppmmax", "rtmed", "rtmin", "rtmax",
                            "npeaks", classnames)

    groupmat <- groupmat[seq(length = num),]
    groupindex <- groupindex[seq(length = num)]

    # Remove groups that overlap with more "well-behaved" groups
    numsamp <- rowSums(groupmat[,(match("npeaks", colnames(groupmat))+1):ncol(groupmat),drop=FALSE])
    uorder <- order(-numsamp, groupmat[,"npeaks"])
    uindex <- xcms:::rectUnique(groupmat[,c("ppmmin","ppmmax","rtmin","rtmax"),drop=FALSE],
                         uorder)
	
	res = list()
	res$groups <- groupmat[uindex,]
	res$groupidx<- groupindex[uindex]
	res$peaks <- peakmatrix
	set.groups.metaboanalyst(res, samp.names, samp.classes)
}

set.groups.metaboanalyst<-function(peaks.result, samp.names, samp.classes) {
    groupmat <- peaks.result$groups;
    groupindex <- peaks.result$groupidx;

    sampnum <- seq(length = length(samp.names))
    intcol <- match("int", colnames(peaks.result$peaks))
    sampcol <- match("sample", colnames(peaks.result$peaks))

    # row is peak, col is sample
    values <- matrix(nrow = length(groupindex), ncol = length(sampnum))

    for (i in seq(along = groupindex)) {
       # for each group, replace multiple peaks from the same sample by their sum
       for(m in sampnum){
            samp.inx<-which(peaks.result$peaks[groupindex[[i]], sampcol]==m)
            if(length(samp.inx)>0){
                 values[i, m] <- sum(peaks.result$peaks[groupindex[[i]][samp.inx], intcol]);
            }else{
                 values[i, m] <- NA;
            }
        }
    }
	values = data.frame(values)
    colnames(values) <- samp.names;

    rownames(values) <- paste(round(groupmat[,paste("ppm", "med", sep="")],5));

	values
}

# auxiliary functions for grouping using metaboanalyst algorithm 

# create matrix as used by MetaboAnalyst
create.metaboanalyst.mat <- function(sample.list){
  mat = matrix();
  allmat = NULL;
  if (ncol(sample.list[[1]]) == 2){
	for (i in 1:length(sample.list)){
		mat = cbind(sample.list[[i]][,1],1000, sample.list[[i]][,2],i)
		allmat = rbind(allmat,mat);
	}
  } else {
	for (i in 1:length(sample.list)){
		mat = cbind(sample.list[[i]][,1],sample.list[[i]][,2], sample.list[[i]][,3],i)
		allmat = rbind(allmat,mat);
	}
  }
  colnames(allmat) = c("ppm","rt","int", "sample")
  allmat  
}


findEqualGreaterM = function(x, values) {
  
  if (!is.double(x)) x <- as.double(x)
  if (!is.double(values)) values <- as.double(values)
  .C("FindEqualGreaterM",
     x,
     length(x),
     values,
     length(values),
     index = integer(length(values)),
     DUP = FALSE, PACKAGE = "xcms")$index + 1
}

rectUnique = function(m, order = seq(length = nrow(m)), xdiff = 0, ydiff = 0) {
  
  nr <- nrow(m)
  nc <- ncol(m)
  if (!is.double(m))
    m <- as.double(m)
  .C("RectUnique",
     m,
     as.integer(order-1),
     nr,
     nc,
     as.double(xdiff),
     as.double(ydiff),
     logical(nrow(m)),
     DUP = FALSE, PACKAGE = "xcms")[[7]]
}
