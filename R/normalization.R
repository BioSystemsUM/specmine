# SAMPLE NORMALIZATION

# method: "sum", "median", "ref.sample", "ref.feature"
# ref: if "ref.sample" or "ref.feature", the reference
# if ref.sample and length(ref) > 1 - the mean of the samples is calculated
# constant: if "sum" method, the constant value

"normalize" = function(dataset, method, ref = NULL, constant = 1000) {
	dataset$data = normalize.samples(dataset$data, method, ref, constant = constant)
	add.desc = paste("Normalization with method", method, sep=" ")
	dataset$description = paste(dataset$description, add.desc, sep="; ")
	dataset
}

"normalize.samples" = function(datamat, method, ref = NULL, constant = 1000) {
  
	if (method == "sum") {
		res = apply(datamat, 2, normalization.sum, constant)			
	} 
  else if (method == "median") {
		res = apply(datamat, 2, normalization.median)
	} 
  else if (method == "ref.sample") {
		if (!is.null(ref) & length(ref) > 0) {
      if (length(ref) == 1)
      {
			  ref.sample = datamat[,ref]
			  res = apply(datamat, 2, normalization.ref.sample, ref.sample)
      }
      else {
        ref.mean = apply(datamat[,ref,drop=F], 1, mean)
        res = apply(datamat, 2, normalization.ref.sample, ref.mean)
      }
		}
		else stop("Reference not defined")
	} 
  else if (method == "ref.feature") {
		if (!is.null(ref) & length(ref) > 0) {
      if (length(ref) == 1) {
        ref.index = which(rownames(datamat) == ref)
			  res = apply(datamat, 2, normalization.ref.feature, ref.index, constant)
      }
      else stop("Reference variable needs to have length 1")
		}
    else stop("Reference not defined")
	}
	else stop("Method not defined")
  
# 	if(method=="ref.feature" && !is.null(ref)) {
# 		index.to.remove = which(rownames(datamat) %in% ref)
#     res = datamat[-index.to.remove,]
#   }
	res
}

"normalization.sum" = function(x, constant = 1000) {
	constant*x/sum(x, na.rm=T);
}

"normalization.median" = function(x) {
	x/median(x, na.rm=T);
}

"normalization.ref.sample" = function(x, ref.sample) {
	x/median(x/ref.sample)
}

"normalization.ref.feature" = function(x, ref.feature, constant = 1000){
	constant*x/x[ref.feature]
}

  
# DATA TRANSFORMATION

# transform data
# experiment: dataset and metadata list
# method: "log" method
"transform.data" = function(dataset, method ="log"){
	if (method == "log") dataset$data = log.transform(dataset$data)
  else if (method == "cubicroot") dataset$data = cubic.root.transform(dataset$data)
  else stop("Method for data transformation not defined")
  
	add.desc = paste("Data transformation with method", method, sep=" ")
	dataset$description = paste(dataset$description, add.desc, sep="; ")
	dataset
}

# generalized log - takem from MetaboAnalyst
"log.transform" = function(datamat) {
  min.val = min(abs(datamat[datamat!=0]))/10
  res = apply(datamat, 2, log.norm, min.val)
}

"cubic.root.transform" = function(datamat) {
  res = apply(datamat, 2, cubic.root)
}

# generalized log, tolerant to 0 and negative values - MetaboAnalyst
"log.norm" = function(x, min.val)
{
  log2((x + sqrt(x^2 + min.val^2))/2)
}

"cubic.root" = function(x) {
  tmp = abs(x)^(1/3)
  tmp[x<0] = -tmp[x<0]
  tmp
}

# SCALING

# scaling dataset
# experiment: dataset and metadata list
# method: "auto", "range", "pareto" methods
"scaling" = function(dataset, method = "auto"){
	dataset$data = scaling.samples(dataset$data, method)
	add.desc = paste("Scaling with method", method, sep=" ")
	dataset$description = paste(dataset$description, add.desc, sep="; ")
	dataset
}


"scaling.samples" = function(datamat, method="auto")
{
  res = datamat
  #for(i in 1:nrow(res))
    if(method == "auto") res = t(apply(res, 1, auto.scale)) #res[i,] = auto.scale(res[i,])
    else if (method == "range") res = t(apply(res, 1, range.scale)) #res[i,] = range.scale(res[i,])
    else if (method == "pareto") res = t(apply(res, 1, pareto.scale)) #res[i,] = pareto.scale(res[i,])
    else if (method == "tointerval") res = t(apply(res, 1, scale.to.interval)) #res[i,] = scale.to.interval(res[i,]) 
    else res[i,] = res[i,] # do nothing
  colnames(res) = colnames(datamat)
  rownames(res) = rownames(datamat)
  res
}

# auto scaling: subtract by mean and divide by sd 
"auto.scale" = function(x){
  (x - mean(x))/sd(x, na.rm=T)
}

"pareto.scale"<-function(x){
  (x - mean(x))/sqrt(sd(x, na.rm=T))
}

"range.scale"<-function(x){
  if(max(x) == min(x)){
    x;
  }
  else{
    (x - mean(x))/(max(x)-min(x));
  }
}

scale.to.interval = function(x, interval = c(0,1)) {
  if(max(x) == min(x)){
    res = x
  }
  else{
    res = (x - min(x))/(max(x)-min(x))
  }

  res = res* (interval[2]-interval[1]) + interval[1]
  res
}