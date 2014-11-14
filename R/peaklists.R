## Functions to read, create and handle peak lists

# Structure: list with one field per sample
# each sample is a data frame representing a peak list (with two columns)

# returns a dataset in standard format from all peaks in all samples
"dataset.from.peaks" = function(sample.list, metadata = NULL, description = "", type = "nmr-peaks") {
  
  merged.peaks = merge.eq.peaks.samplelist(sample.list)
  samples.df = get.all.intensities(merged.peaks)
  create.dataset(as.matrix(samples.df), metadata = metadata, type = type, 
                 description = description)

}

# READING

# reads csv files, each with a sample; 
# filenames - list of file names of the files to read
# returns list of data frames
"read.multiple.csvs" = function(filenames, ext = ".csv", ...)
{
  sampleList = list()
  sampleNames = c()
  snames <- gsub("\\.[^.]*$", "", basename(filenames));
  for (i in 1:length(filenames)) {
    print(paste("Reading sample ", filenames[i]))
    sampleList[[i]] = read.csv(paste(filenames[i], ext, sep=""), ...)
  }
  sampleNames = snames
  names(sampleList) = sampleNames
  sampleList
}

# reads list of CSV files from a given folder
# returns list of data frames
"read.csvs.folder" = function(foldername, ...)
{
  files<-dir(foldername, pattern=".[Cc][Ss][Vv]$", recursive=T, full.name=TRUE)
  sampleList = read.multiple.csvs(files, ext= "", ...);
  sampleList
}

# GETTING INFO

# counts number of peaks in a sample (given its index)
"peaks.per.sample" = function(sample.list, sample.index)
{
  nrow(sample.list[[sample.index]])
}

# counts the number of peaks in each sample in the peaklist
"peaks.per.samples" = function(sample.list)
{
  res = c()
  for(i in 1:length(sample.list)) res[i] = peaks.per.sample(sample.list, i)
  res
}

# finds samples that have the same peak values- x and y (equal data frames)
"find.equal.samples" = function(sample.list)
{
  require(compare)
  eq1 = c()
  eq2 = c()
  for (i in 1:(length(sample.list)-1))
  {
    for(j in (i+1):length(sample.list)) 
    {
      res = compare(sample.list[[i]], sample.list[[j]])
      if(res$result == T) 
      {
        eq1 = c(eq1, names(sample.list)[i])
        eq2 = c(eq2, names(sample.list)[j])
      }
    }
  }
  data.frame(eq1, eq2)
}

# get the full list of frequencies from all samples in a list
"get.overall.freq.list" = function(sample.list)
{
  res = sample.list[[1]][[1]]
  for(i in 2:length(sample.list))
    res = union(res, sample.list[[i]][[1]])
  sort(res)
}

# get the value of an internsity given the sample (data frame) and the frequency
"get.intensity" = function(sample.df, freq, tolerance = 0.001)
{
  cond = (sample.df[[1]] > freq - tolerance) & (sample.df[[1]] < freq + tolerance)
  if (any(cond)) res = sample.df[cond,][[2]]
  else res = NA
  res
}

# gets all intensities for a sample list; result is a data frame
"get.all.intensities" = function(sample.list, tol = 0.001)
{
  all.freqs = get.overall.freq.list(sample.list)
  
  for(i in 1:length(sample.list))
  {
    intens.vals = c()
    for(k in 1:length(all.freqs)) 
    {
      intens.vals[k] = get.intensity(sample.list[[i]], all.freqs[k], tol)
    }
    if (i==1) res.df = data.frame(intens.vals)
    else res.df = cbind(res.df, intens.vals)
  }
  names(res.df) = names(sample.list)
  rownames(res.df) = all.freqs
  res.df
}

"remove.peaks.interval" = function(sample.df, peak.val.min, peak.val.max) {
  subset(sample.df, ppm < peak.val.min | ppm > peak.val.max)
}

"remove.peaks.interval.sample.list" = function(sample.list, peak.val.min, peak.val.max) {
  sample.list.res = list()
  for(i in 1:length(sample.list)) {
    sample.list.res[[i]] = remove.peaks.interval(sample.list[[i]], peak.val.min, peak.val.max)
  }
  names(sample.list.res) = names(sample.list)
  sample.list.res
}

# functions working over the list of all intensities
"get.peak.values" = function(samples.df, peak.val)
{
  index = which(rownames(samples.df) == peak.val)
  values = c()
  for (s in (samples.df[index,])) {
    values = c(values, s)
  }
  values
}

"values.per.peak" = function(samples.df)
{
  res = c()
  for(i in 1:nrow(samples.df)) res[i] = sum(!is.na(samples.df[i,]))
  res
}

"values.per.sample" = function(samples.df)
{
  res = c()
  for(i in 1:ncol(samples.df)) res[i] = sum(!is.na(samples.df[,i]))
  res
}

# PROCESSING

# merge peaks with equal frequencies (ppm) in a sample given by a data frame
# intensity values are summed
"merge.equal.peaks" = function(sample.df, tolerance = 0.0)
{
  d = diff(sample.df[[1]])
  indexes = which(d <= tolerance)
  if (length(indexes) != 0){
	new.sample.df = sample.df[-(indexes+1),]
	new.sample.df[[2]] = sum.vec(sample.df[[2]], indexes)
  }
  else {
	new.sample.df = sample.df
  }

  new.sample.df
}

"sum.vec" = function(orig.vec, indexes)
{
  newvec = c()
  for (i in 1:length(orig.vec)) newvec[i] = orig.vec[i]
  
  for(k in length(indexes):1)
  {
    newvec[indexes[k]] = newvec[indexes[k]] + newvec[indexes[k]+1]
    newvec = newvec[-(indexes[k]+1)]
  }
  newvec
}

# merge peaks with equal freqs (ppm) in all samples of a list 
"merge.eq.peaks.samplelist" = function(sample.list, tolerance = 0.0)
{
  newlist = list()
  for(i in 1:length(sample.list))
  {
    newlist[[i]] = merge.equal.peaks(sample.list[[i]], tolerance)
  }
  names(newlist) = names(sample.list)
  newlist
}

# create a matrix with all values for all samples (as used by MetaboAnalyst)
"create.full.matrix" = function(sample.list)
{
  resmat = NULL
  for(i in 1:length(sample.list))
  {
    for(j in 1:nrow(sample.list[[i]]))
      resmat = rbind(resmat, cbind(sample.list[[i]][j,], i))
  }
  colnames(resmat) = c("ppm","int","sample")
  resmat
}

