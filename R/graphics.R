
### BOXPLOT of each variable ####

"boxplot.variables" = function(dataset, variables = NULL, samples = NULL, horizontal = T, 
                               col = "lightblue", nchar.label = 10, cex.axis = 0.8, ...) {
  
  if (is.null(variables)) { # assume all
    variables = rownames(dataset$data)
  } 
  
  if (is.null(samples)) {
    samples = colnames(dataset$data)
  }
  if (is.numeric(variables))
    names.short = substr(get.x.values.as.text(dataset)[variables], 1, nchar.label)
  else 
    names.short = substr(variables, 1, nchar.label)

  boxplot(t(dataset$data[variables,samples]), names = names.short, 
          horizontal = horizontal, las = 2, 
          col = col, cex.axis = cex.axis, ...)
}


##############################SPECTRA PLOTS###############################

# samples - list of samples to plot
# variable.bounds - interval of x values to plot: [1] - minimum value; [2] - maximum value
# lty, lwd, col - parameters to pass to matplot
# ... - extra parameters passed to matplot function
"plot.spectra.simple" = function(dataset, samples = NULL, variable.bounds = NULL, xlab = NULL,
                               ylab = NULL, lty = 1, lwd = 1, col = 1, reverse.x = F, ...) {
  
  if (is.null(xlab)) xlab = get.x.label(dataset)
  if (is.null(ylab)) ylab = get.value.label(dataset)
  
  if (!is.null(dataset$labels$x) && dataset$labels$x == "mz/rt"){
	if (is.null(variable.bounds)){
		variables = 1:length(get.x.values.as.text(dataset))
		vars = variables
	} else {
		x.vars = as.numeric(gsub("/.*", '', get.x.values.as.text(dataset)))
		variables = which(x.vars > variable.bounds[1] & x.vars < variable.bounds[2])
		vars = x.vars[variables]
	}
  } else {
	  if (is.null(variable.bounds)){
		variables = rownames(dataset$data)
		vars = variables
	  } 
	  else {
		x.vars = get.x.values.as.num(dataset)
		variables = rownames(dataset$data)[x.vars > variable.bounds[1] & x.vars < variable.bounds[2]] 
		vars = variables
	  }
  }
  if (reverse.x) xlim = c( max(as.numeric(vars)), min(as.numeric(vars)) )
  else xlim = range(as.numeric(vars))
  
  if (is.null(samples)){
    samples = colnames(dataset$data)
  } 

  matplot(vars, dataset$data[variables,samples,drop=F], type="l", lty=lty, col = col,
            xlab = xlab, ylab = ylab, xlim = xlim, ...)
}


# plots spectra with colors per groups, given by a metadata variable

# legend.place = "none" if no legend or other accepted string in legend function
# 
"plot.spectra" = function(dataset, column.class, func = NULL, samples = NULL, 
                        variable.bounds = NULL, xlab = NULL, ylab = NULL, lty = 1,
                        legend.place = "topright", cex = 0.8, reverse.x = F, ...) {
  
  if (is.null(xlab)) xlab = get.x.label(dataset)
  if (is.null(ylab)) ylab = get.value.label(dataset)
  
 if (!is.null(dataset$labels$x) && dataset$labels$x == "mz/rt"){
	if (is.null(variable.bounds)){
		variables = 1:length(get.x.values.as.text(dataset))
		vars = variables
	} else {
		x.vars = as.numeric(gsub("/.*", '', get.x.values.as.text(dataset)))
		variables = which(x.vars > variable.bounds[1] & x.vars < variable.bounds[2])
		vars = x.vars[variables]
	}
  } else {
	  if (is.null(variable.bounds)){
		variables = rownames(dataset$data)
		vars = variables
	  } 
	  else {
		x.vars = get.x.values.as.num(dataset)
		variables = rownames(dataset$data)[x.vars > variable.bounds[1] & x.vars < variable.bounds[2]] 
		vars = variables
	  }
  }
  
  if (reverse.x) xlim = c( max(as.numeric(vars)), min(as.numeric(vars)) )
  else xlim = range(as.numeric(vars))
	
  if (is.null(samples)){
		samples = colnames(dataset$data)
		metadata = dataset$metadata[,column.class]
	} 
  else {
		metadata = factor(dataset$metadata[samples, column.class])
	}
	if (is.null(func)){
		matplot(vars, dataset$data[variables,samples], type="l", lty=lty, col=as.integer(metadata), 
            xlab = xlab, ylab = ylab, xlim = xlim, ...)
    if (legend.place != "none")
		  legend(legend.place, levels(metadata), cex=cex, fill = sort(as.integer(factor(levels(metadata))))) 
	} 
  else {
		aggregate.result = aggregate(t(dataset$data[variables,samples]), by = list(metadata), func)
		matplot(vars, t(aggregate.result[-1]), type = "l", lty=1, col=as.integer(aggregate.result[,1]), 
            xlab = xlab, ylab = ylab, xlim = xlim, ...)
		if (legend.place != "none")
		  legend(legend.place, levels(metadata), cex=cex, fill = sort(as.integer(factor(levels(metadata)))))
	}
}

