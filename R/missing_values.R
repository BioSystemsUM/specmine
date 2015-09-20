# diagnostics

"count.missing.values" = function(dataset)
{
  sum(is.na(dataset$data))
}

"count.missing.values.per.sample" = function(dataset, remove.zero = T) {
  res = apply(dataset$data, 2, function(x) sum(is.na(x)))
  if (remove.zero) res[res > 0]
  else res
}

"count.missing.values.per.variable" = function(dataset, remove.zero = T) {
  "count.na" = function(x) sum(is.na(x))
  res = apply(dataset$data, 1, count.na)
  if (remove.zero) res[res > 0]
  else res
}


# missing values imputation

# method: "value", "mean", "median" or "knn" imputation method
# value: if "value" method selected, the value that will replace NAs
# k: if "knn" method selected, the number of neighbors
"missingvalues.imputation" = function(dataset, method = "value", value = 0.0005, k = 5){
	if (method == "value"){
		dataset = impute.nas.value(dataset, value)
	} 
  else if (method == "mean"){
		dataset = impute.nas.mean(dataset)
	} 
  else if (method == "median"){
		dataset = impute.nas.median(dataset)
	} 
  else if (method == "knn"){
		dataset = impute.nas.knn(dataset, k)
	} 
  else if (method == "linapprox"){
		dataset = impute.nas.linapprox(dataset)
	}
  add.desc = paste("Missing value imputation with method", method, sep=" ")
	dataset$description = paste(dataset$description, add.desc, sep="; ")
  dataset
}


impute.nas.linapprox = function(dataset){
  hyper.object = convert.to.hyperspec(dataset)
	linapprox.res = spc.NA.linapprox(hyper.object)
  res.dataset = convert.from.hyperspec(linapprox.res)
  res.dataset$type = dataset$type
  res.dataset
}

"impute.nas.value" = function(dataset, value)
{
  dataset$data[is.na(dataset$data)] = value
  dataset
}

#taken from metaboanalyst
"impute.nas.mean" = function(dataset){
  temp = apply(dataset$data, 1, function(x){
					if(sum(is.na(x))>0){
						x[is.na(x)] = mean(x, na.rm=T);
					}
					x;
				})
  dataset$data = t(temp)
  dataset
}

#taken from metaboanalyst
"impute.nas.median" = function(dataset) {
  temp = apply(dataset$data, 1, function(x){
					if(sum(is.na(x))>0){
						x[is.na(x)] = median(x,na.rm=T);
					}
					x;
				})
  dataset$data = t(temp)
  dataset
}

#taken from metaboanalyst - uses impute package from bioconductor

"impute.nas.knn" = function(dataset, k = 10, ...){
  dataset$data = impute::impute.knn(dataset$data, ...)
  dataset
}

