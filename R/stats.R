# applies a function to the values of each variable
# fn.to.apply - function to apply (e.g. mean, max, min)
# variables - allows to define which variables to calculate the stats (if numbers, indexes are assumed)
# variable.bounds - allow to define an interval of variables (if numeric)
# samples - if defined restricts the application to a given set of samples 
"apply.by.variable" = function(dataset, fn.to.apply, variables = NULL, variable.bounds = NULL, 
                               samples = NULL, ...) {
  
  if (is.null(variables)) {
    if (is.null(variable.bounds)){
      variables = rownames(dataset$data)
    } 
    else {
      x.vars = get.x.values.as.num(dataset)
      variables = rownames(dataset$data)[x.vars > variable.bounds[1] & x.vars < variable.bounds[2]] 
    }  
  }  
  if (is.null(samples)) {
    samples = colnames(dataset$data)
  }
  
  apply(dataset$data[variables,samples,drop=F], 1, fn.to.apply, ...)
}

"apply.by.sample" = function(dataset, fn.to.apply, samples = NULL, ...) {
  if (is.null(samples)) {
    samples = colnames(dataset$data)
  }
  apply(dataset$data[,samples,drop=F], 2, fn.to.apply, ...)
}

"stats.by.variable" = function(dataset, variables = NULL, variable.bounds = NULL) {
  apply.by.variable(dataset, summary, variables, variable.bounds)
}

"stats.by.sample" = function(dataset, samples = NULL) {
  apply.by.sample(dataset, summary, samples)
}

"apply.by.group" = function(dataset, fn.to.apply, metadata.var, var.value) {
  indexes = which(dataset$metadata[,metadata.var] %in% var.value)
  apply.by.sample(dataset, fn.to.apply, indexes)
}