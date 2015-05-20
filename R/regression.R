"linregression.onevar" = function(dataset, x.val, metadata.vars, combination) {
  values = get.data.values(dataset, x.val)
  sub.df = dataset$metadata[,metadata.vars]
  sub.df = cbind(values,sub.df)
  terms = names(sub.df)[2:ncol(sub.df)]
  terms = cbind(terms, combination)
  reslm = lm(reformulate(terms, "values"), data = sub.df)
  lm.summary = summary(reslm)
  lm.summary
}

"linreg.all.vars" = function(dataset, metadata.vars, combination)
{
  m = vector("list",nrow(dataset$data))
  for(i in 1:nrow(dataset$data))
    m[[i]] = linregression.onevar(dataset, rownames(dataset$data)[i], metadata.vars, combination)
  names(m) = rownames(dataset$data)  
  m
}

linreg.coef.table = function(linreg.results, write.file = F, file.out = "linreg-coefs.csv"){
  num_vars = dim(linreg.results[[1]]$coefficients)[1]
  m = matrix(NA, length(linreg.results), num_vars)
  rownames(m) = names(linreg.results)
  for (i in 1:length(linreg.results)){
    m[i,] = linreg.results[[i]]$coefficients[,1]
  }
  coef.table = as.data.frame(m)
  colnames(coef.table) = trim(rownames(linreg.results[[1]]$coefficients))
  if (write.file) write.csv(coef.table, file = file.out)
  coef.table  
}

linreg.pvalue.table = function(linreg.results, write.file = F, file.out = "linreg-pvalues.csv"){
  num_vars = dim(linreg.results[[1]]$coefficients)[1]
  m = matrix(NA, length(linreg.results), num_vars)
  rownames(m) = names(linreg.results)
  for (i in 1:length(linreg.results)){
    m[i,] = linreg.results[[i]]$coefficients[,4]
  }
  pv.table = as.data.frame(m)
  colnames(pv.table) = trim(rownames(linreg.results[[1]]$coefficients))
  if (write.file) write.csv(coef.table, file = file.out)
  pv.table  
}

linreg.rsquared = function(linreg.results, write.file = F, file.out = "linreg-rsquared.csv"){
  m = matrix(NA, length(linreg.results), 2)
  rownames(m) = names(linreg.results)
  for (i in 1:length(linreg.results)){
    m[i,1] = linreg.results[[i]]$r.sq
    m[i,2] = linreg.results[[i]]$adj.r.sq
  }
  rsq.table = as.data.frame(m)
  colnames(rsq.table) = c("r.squared", "adj.r.squared")
  if (write.file) write.csv(rsq.table, file = file.out)
  rsq.table  
}