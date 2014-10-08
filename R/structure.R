## FUNCTIONS TO DEFINE AND QUERY DATASET STRUCTURE

list.of.spectral.types = c("nmr-spectra","ir-spectra", "uvv-spectra", "raman-spectra", "fluor-spectra")

list.of.allowed.types = c(list.of.spectral.types, "nmr-peaks", "concentrations", "undefined")


# function to create a dataset from existing objects 
# datamatrix - matrix with numerical data; rows are assumed to be variables and columns assumed to be samples
# type - type of data: can be one of the following: "nmr-spectra", "nmr-peaks", "ir-spectra", "uvv-spectra", 
# "concentrations", "undefined", ...

"create.dataset" = function(datamatrix, type = "undefined", metadata = NULL, description = "", 
                            sample.names = NULL, x.axis.values = NULL, 
                            label.x = NULL, label.values = NULL) {
  
  if (is.null(datamatrix))
    stop("Invalid argument: datamatrix is null")
  
  if (!is.matrix(datamatrix)) {
    if (is.data.frame(datamatrix)) {
      warning("datamatrix is data.frame; converting to matrix")
      datamatrix = as.matrix(datamatrix)
    }
    else stop("Invalid argument: datamatrix is not matrix or data.frame")
  }
  if (!is.numeric(datamatrix))
    stop("datamatrix is not numeric")
  
  if (! type %in% list.of.allowed.types) 
    stop("Type of data is not allowed")
  
  if (!is.null(metadata)) {
    if (nrow(metadata) != ncol(datamatrix))
      stop("Number of columns in data matrix not the same as number of rows in metadata")
    if (!is.data.frame(metadata)){
      if (is.matrix(metadata)) metadata = as.data.frame(metadata)
      else stop("metadata is not matrix or data.frame")
    }
  }
  else warning("Metadata is null; dataset will still be created with empty metadata")
  
  if (!is.null(label.x) | !is.null(label.values) ) 
    labels = list(x = label.x, val = label.values)
  else {
    labels = NULL
    warning("Labels are null")
  }
  
  if (!is.null(sample.names)) {
    if (length(sample.names) != ncol(datamatrix) )
      stop("Number of columns in data matrix not the same as length of sample names vector")
    colnames(datamatrix) = sample.names

    if (!is.null(metadata)) {
      if(length(sample.names) != nrow(metadata))
        stop("Number of rows in metadata not the same as length of sample names vector")
      rownames(metadata) = sample.names
    }
  }
  else {
    if (is.null(colnames(datamatrix))) { # default names will be row/col numbers
      warning("Sample names not specified; will be assumed as sequential numbers")
      colnames(datamatrix) = as.character(1:ncol(datamatrix))
      if (!is.null(metadata)) rownames(metadata) = as.character(1:nrow(metadata))
    }
  }
  
  if (!is.null(x.axis.values)) {
    if (length(x.axis.values) != nrow(datamatrix))
      stop("Number of rows in data matrix not the same as length of x axis values vector")
    if (type %in% list.of.spectral.types & any(is.na(as.numeric(x.axis.values))) )
      stop("Invalid non numeric values for variable names in x.axis.values parameter (given spectral type)")
    rownames(datamatrix) = as.character(x.axis.values)
  }
  else {
    if (is.null(rownames(datamatrix))) {
      warning("data variable names not specified; will be assumed as sequential numbers")
      rownames(datamatrix) = as.character(1:nrow(datamatrix))
    }
    else
      if (type %in% list.of.spectral.types) 
        if(any(is.na(as.numeric(rownames(datamatrix)))) )
          stop("Invalid non numeric values for variable names in rownames of matrix (given spectral type)")
  }
  
  dataset = list(data = datamatrix, type = type, description = description, metadata = metadata, labels = labels)
  # make sure sample names are the same in data and metadata
  if (!is.null(dataset$metadata)) rownames(dataset$metadata) = colnames(dataset$data)
  
  dataset
}

"check.dataset" = function(dataset)
{
  if (is.null(dataset$data)) 
    stop("Invalid dataset: Data matrix is null")
  
  if (!is.null(dataset$metadata)) {
    if (nrow(dataset$metadata) != ncol(dataset$data) )
      stop("Invalid dataset: Number of columns in data matrix not the same as number of rows in metadata")
  }
  else warning("Metadata is null")
  
  if (!dataset$type %in% list.of.allowed.types) stop("Type of data is not allowed")
  
  if (dataset$type %in% list.of.spectral.types) 
    if (any(is.na(as.numeric(rownames(dataset$data)))) )
      stop("Invalid non numeric values for variable names in rownames of matrix (given spectral type)")
  
  cat("Valid dataset\n")
  res = TRUE
}

# provides a summary of the dataset, printing its main features
# stats - if TRUE prints some global statistics of the data values
"sum.dataset" = function(dataset, stats = T)
{
  cat("Dataset summary:\n")
  check.dataset(dataset)
  cat ("Description: ", dataset$description, "\n")
  cat("Type of data: ", dataset$type, "\n")
  cat("Number of samples: ", ncol(dataset$data), "\n")
  cat("Number of data points", nrow(dataset$data), "\n")
  if (!is.null(dataset$metadata))
    cat("Number of metadata variables: ", ncol(dataset$metadata), "\n")
  if (!is.null(dataset$labels)) {
    if (!is.null(dataset$labels$x)) 
      cat("Label of x-axis values: ", as.character(dataset$labels$x), "\n")
    if (!is.null(dataset$labels$val)) 
      cat("Label of data points: ", as.character(dataset$labels$val), "\n")
  }
  if (stats) {
    cat("Number of missing values in data: ", sum(is.na(dataset$data)), "\n")
    cat("Mean of data values: ", mean(dataset$data, na.rm= T), "\n")
    cat("Median of data values: ", median(dataset$data, na.rm = T), "\n")
    cat("Standard deviation: ", sd(dataset$data, na.rm = T), "\n")
    cat("Range of values: ", range(dataset$data, na.rm = T), "\n")
    cat("Quantiles:", "\n")
    print(quantile(dataset$data, na.rm=T))
  }
}

# QUERY functions
# functions to access data from a dataset

# returns data matrix
"get.data" = function(dataset)
{
  dataset$data
}

# returns data matrix as a data frame
"get.data.as.df" = function(dataset)
{
  as.data.frame(dataset$data)
}

# returns metadata (data frame)
"get.metadata" = function(dataset)
{
  dataset$metadata
}

# returns values of a metadata variable
# var - index or name of the metadata variable
"get.metadata.var" = function(dataset, var)
{
  dataset$metadata[,var]
}

"num.samples" = function(dataset)
{
  ncol(dataset$data)
}

"get.sample.names" = function(dataset)
{
  sample.names = colnames(dataset$data)
  sample.names
}

"num.x.values" = function(dataset)
{
  nrow(dataset$data)
}

"get.x.values.as.text" = function(dataset)
{
  x.values = rownames(dataset$data)
  as.character(x.values)
}

"get.x.values.as.num" = function(dataset)
{
  x.values = rownames(dataset$data)
  res = as.numeric(x.values)
  if (any(is.na(res))) stop("Variable labels are not all numeric")
  res
}

"get.x.label" = function(dataset) {
  if (is.null(dataset$labels) | is.null(dataset$labels$x)) return ("")
  else return (dataset$labels$x)
}


"get.value.label" = function(dataset) {
  if (is.null(dataset$labels) | is.null(dataset$labels$val)) return ("")
  else return (dataset$labels$val)
}

"get.type" = function(dataset) {
  dataset$type
}

# specifies if a dataset is from spectral data where x.values are numeric
"is.spectra" = function(dataset) {
  dataset$type %in% list.of.spectral.types
}

# returns a data value given the x axis labes (as index or name) and the sample (as index or name)
"get.data.value" = function(dataset, x.axis.val, sample, by.index = F) {
  if (!by.index) {
    x.axis.val = as.character(x.axis.val)
    x.axis.index = which(rownames(dataset$data) == x.axis.val)
  }
  else {
    x.axis.index = x.axis.val
  }
  dataset$data[x.axis.index, sample]
}

# can use both indexes or names
"get.metadata.value" = function(dataset, variable, sample)
{
  dataset$metadata[sample, variable]
}

# returns values of all samples given a set of x axis names (or indexes of by,index is T)
"get.data.values" = function(dataset, x.axis.val, by.index = FALSE)
{
  if (!by.index) {
    if (length(x.axis.val) >= 1) {
      x.axis.val = as.character(x.axis.val)
      x.axis.indexes = which(rownames(dataset$data) %in% x.axis.val)
    }
    else 
      stop("Incorrect parameter x.axis.val: length not >= 1")
  }
  else {
    x.axis.indexes = x.axis.val
  }
  
  dataset$data[x.axis.indexes,]
}

# returns indexes corresponding to a vector of x-values (assuming numerical values - spectra)
"x.values.to.indexes" = function(dataset, x.values)
{
  x.values.ds = get.x.values.as.num(dataset)
  indexes = which(x.values.ds %in% x.values)
  indexes
}

# returns indexes corresponding to an interval of x-values (assuming numerical values - spectra)
"xvalue.interval.to.indexes" = function(dataset, min.value, max.value) {
  x.values = get.x.values.as.num(dataset)
  indexes = which(x.values >= min.value & x.values <= max.value)
  indexes
}

"indexes.to.xvalue.interval" = function(dataset, indexes) {
  x.values = get.x.values.as.num(dataset)
  x.val.inds = x.values[indexes]
  c(min(x.val.inds), max(x.val.inds))
}


# UPDATE functions

"set.metadata" = function(dataset, new.metadata)
{
  if (nrow(new.metadata) != ncol(dataset$data)) 
    stop("Number of columns in data matrix not the same as number of rows in metadata")
  if (!is.data.frame(new.metadata)){
      if (is.matrix(new.metadata)) new.metadata = as.data.frame(new.metadata)
      else stop("metadata is not matrix or data.frame")
  }
  
  dataset$metadata = new.metadata
  dataset
}

"set.x.values" = function(dataset, new.x.values, new.x.label = NULL)
{
  if (length(new.x.values) != nrow(dataset$data) )
    stop("Length of new vector is not consistent with dataset")
  rownames(dataset$data) = as.character(new.x.values)
  if (!is.null(new.x.label))
    dataset = set.x.label(dataset, new.x.label)
  dataset
}

"set.x.label" = function(dataset, new.x.label)
{
  if (!is.null(dataset$label))
    dataset$labels$x = new.x.label 
  else {
    dataset$labels = list()
    dataset$labels$x = new.x.label
  }
  dataset
}

"set.value.label" = function(dataset, new.val.label)
{
  if (!is.null(dataset$label))
    dataset$labels$val = new.val.label 
  else {
    dataset$labels = list()
    dataset$labels$val = new.val.label
  }
  dataset
}

"set.sample.names" = function(dataset, new.sample.names)
{
  if (length(new.sample.names) != ncol(dataset$data))
    stop("Length of new sample names not consistent with dataset dimensions")
  
  colnames(dataset$data) = new.sample.names
  rownames(dataset$metadata) = new.sample.names
  
  dataset
}

"replace.data.value" = function(dataset, x.axis.val, sample, new.value, by.index = F) {
  if (!by.index) {
    x.axis.val = as.character(x.axis.val)
    x.axis.index = which(rownames(dataset$data) == x.axis.val)
  }
  else {
    x.axis.index = x.axis.val
  }
  dataset$data[x.axis.index, sample] = new.value
  dataset
}

"replace.metadata.value" = function(dataset, variable, sample, new.value)
{
  dataset$metadata[sample, variable] = new.value
  dataset
}

"convert.to.factor" = function(dataset, metadata.var)
{
  dataset$metadata[,metadata.var] = factor(dataset$metadata[,metadata.var])
  dataset
}

# MERGE DATASETS
# merges two datasets; data and metadata variables are assumed to be the same and kept from dataset1
# samples from both datasets are merged
"merge.datasets" = function(dataset1, dataset2)
{
  if (ncol(dataset1$metadata) != ncol(dataset2$metadata))
    stop("Different number of metadata variables")
  if (nrow(dataset1$data) != nrow(dataset2$data))
    stop("Different number of data variables")
  
  dataset1$data = cbind(dataset1$data, dataset2$data)
  dataset1$metadata = rbind(dataset1$metadata, dataset2$metadata)
  
  dataset1
}