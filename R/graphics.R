
### BOXPLOT of each variable ####

"boxplot_variables" = function(dataset, variables = NULL, samples = NULL, horizontal = TRUE, 
                               col = "lightblue", nchar.label = 10, cex.axis = 0.8, ...) {
  
  if (is.null(variables)) { # assume all
    variables = rownames(dataset$data)
  } 
  
  if (is.null(samples)) {
    samples = colnames(dataset$data)
  }
  if (is.numeric(variables))
    names.short = substr(get_x_values_as_text(dataset)[variables], 1, nchar.label)
  else 
    names.short = substr(variables, 1, nchar.label)

  if (length(variables) > 1) {
      boxplot(t(dataset$data[variables,samples]), names = names.short, 
          horizontal = horizontal, las = 2, 
          col = col, cex.axis = cex.axis, ...)
  }
  else {
      boxplot(dataset$data[variables,samples], xlab = names.short,  
            horizontal = FALSE, las = 2, 
            col = col, cex.axis = cex.axis, ...)
  }
}

boxplot_vars_factor = function(dataset, meta.var, variables = NULL, samples = NULL, 
                               horizontal = FALSE, nchar.label = 10, col = NULL,
                               vec.par = NULL, cex.axis = 0.8, ylabs = NULL, ...)
{
  if (is.null(variables)) { # assume all
    variables = rownames(dataset$data)
  } 
  
  if (is.null(samples)) {
    samples = colnames(dataset$data)
  }
  
  if (is.numeric(variables))
    names.short = substr(get_x_values_as_text(dataset)[variables], 1, nchar.label)
  else 
    names.short = substr(variables, 1, nchar.label)
  
  if (is.null(vec.par))
    vec.par = c(length(variables), 1)
  
  opar=par(mfrow = vec.par)
  on.exit(par(opar))
  for (i in 1:length(variables)) {
    if (is.null(col)) coli = i+1
    else coli = col
	if (!is.null(ylabs)){
		ylab = ylabs[i]
	} else {
		ylab = NULL
	}
    boxplot(dataset$data[variables[i],samples] ~ dataset$metadata[,meta.var], 
            horizontal = horizontal, las = 2, main = names.short[i], 
            col = coli, cex.axis = cex.axis, ylab = ylab, ...)
  }
  opar=par(mfrow = c(1,1))
  on.exit(par(opar))
}

plotvar_twofactor = function(dataset, variable, meta.var1, meta.var2, colour = "darkblue", title = "", 
                             xlabel = NULL, ylabel = NULL)
{

  df = data.frame(dataset$data[variable,], dataset$metadata[,meta.var1], 
                  dataset$metadata[,meta.var2])
  if (is.numeric(variable)) n1 = rownames(dataset$data)[variable]
  else n1 = variable
  if (is.numeric(meta.var1)) n2 = colnames(dataset$metadata)[meta.var1]
  else n2 = meta.var1
  if (is.numeric(meta.var2)) n3 = colnames(dataset$metadata)[meta.var2]
  else n3 = meta.var2
  colnames(df) = c("name1", "name2", "name3")

  g = ggplot2::ggplot(data = df)
  g = g + ggplot2::aes_string('name2', 'name1')
  g = g + ggplot2::geom_boxplot(fill= colour)
  g = g + ggplot2::facet_grid(. ~ name3)
  if (is.null(xlabel)) g = g + ggplot2::xlab(n2)
  else g = g + ggplot2::xlab(xlabel)
  if (is.null(ylabel)) g = g + ggplot2::ylab(n1)
  else g = g + ggplot2::ylab(ylabel)
  if (title != "") g = g + ggplot2::ggtitle(title)
  g
}

##############################SPECTRA PLOTS###############################

# samples - list of samples to plot
# variable.bounds - interval of x values to plot: [1] - minimum value; [2] - maximum value
# lty, lwd, col - parameters to pass to matplot
# ... - extra parameters passed to matplot function
"plot_spectra_simple" = function(dataset, samples = NULL, variable.bounds = NULL, xlab = NULL,
                               ylab = NULL, lty = 1, lwd = 1, col = 1, reverse.x = FALSE, ...) {
  
  if (is.null(xlab)) xlab = get_x_label(dataset)
  if (is.null(ylab)) ylab = get_value_label(dataset)
  
  if (!is.null(dataset$labels$x) && !is.expression(dataset$labels$x) && dataset$labels$x == "mz/rt"){
	if (is.null(variable.bounds)){
		variables = 1:length(get_x_values_as_text(dataset))
		vars = as.numeric(gsub("/.*", '', get_x_values_as_text(dataset)))

	} else {
		x.vars = as.numeric(gsub("/.*", '', get_x_values_as_text(dataset)))
		variables = which(x.vars > variable.bounds[1] & x.vars < variable.bounds[2])
		vars = x.vars[variables]
	}
  } else {
	  if (is.null(variable.bounds)){
		variables = rownames(dataset$data)
		vars = variables
	  } 
	  else {
		x.vars = get_x_values_as_num(dataset)
		variables = rownames(dataset$data)[x.vars > variable.bounds[1] & x.vars < variable.bounds[2]] 
		vars = variables
	  }
  }
  if (reverse.x) xlim = c( max(as.numeric(vars)), min(as.numeric(vars)) )
  else xlim = range(as.numeric(vars))
  
  if (is.null(samples)){
    samples = colnames(dataset$data)
  } 

  matplot(vars, dataset$data[variables,samples,drop=FALSE], type="l", lty=lty, col = col,
            xlab = xlab, ylab = ylab, xlim = xlim, ...)
}


# plots spectra with colors per groups, given by a metadata variable

# legend.place = "none" if no legend or other accepted string in legend function
# 
"plot_spectra" = function(dataset, column.class, func = NULL, samples = NULL, 
                        variable.bounds = NULL, xlab = NULL, ylab = NULL, lty = 1,
                        legend.place = "topright", cex = 0.8, reverse.x = FALSE, ...) {
  
  if (is.null(xlab)) xlab = get_x_label(dataset)
  if (is.null(ylab)) ylab = get_value_label(dataset)
  
 if (!is.null(dataset$labels$x) && !is.expression(dataset$labels$x) && dataset$labels$x == "mz/rt"){
	if (is.null(variable.bounds)){
		variables = 1:length(get_x_values_as_text(dataset))
		vars = as.numeric(gsub("/.*", '', get_x_values_as_text(dataset)))
	} else {
		x.vars = as.numeric(gsub("/.*", '', get_x_values_as_text(dataset)))
		variables = which(x.vars > variable.bounds[1] & x.vars < variable.bounds[2])
		vars = x.vars[variables]
	}
  } else {
	  if (is.null(variable.bounds)){
		variables = rownames(dataset$data)
		vars = variables
	  } 
	  else {
		x.vars = get_x_values_as_num(dataset)
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

#################################
############2D PLOT##############
#################################

#Function to provide the signal-to-noise ratio of a spectra
snr_spectra <- function(spec){
  threshold <- mean(base::Filter(isPositive,spec))/(sd(base::Filter(isPositive,spec[which(apply(spec,1,cv)>15), which(apply(spec,2,cv)>15)])))
  threshold
}

#Function to provide an ordered data frame of each sample spectra SNR
snr_all <- function(spect, metadata = NULL) {
  tresholds <- unlist(lapply(spect$data,snr_spectra))
  samples <- names(spect$data)
  if (!is.null(metadata)) {
    meta_data <- as.factor(spect$metadata[samples,metadata])
    res <- data.frame(SNR = tresholds, Meta = meta_data, row.names = samples)
    res <- res[order(-res$SNR,res$Meta),]
  } else {
    res <- data.frame(SNR = tresholds, Samples = samples)
    res <- res[order(-res$SNR),]
  }
  res
}

#Function to plot in an interactive way some samples of the specmine_2d_dataset
#if none given, the two higher and two lowest SNR spectra are plotted
plot_2d_spectra <- function(specmine_2d_dataset, title_spectra = "", meta = NULL,spec_samples = NULL) {
  clrs <- c("Greys","YlGnBu","Greens","YlOrRd","Bluered","RdBu","Reds","Blues","Picnic","Rainbow","Portland","Jet","Hot","Blackbody","Earth","Electric","Viridis","Cividis")
  #colorscale2 <- list(c(0, 1), c("tan", "blue"))
  
  if (!is.null(spec_samples)) {
    if (typeof(spec_samples)=="character") {
      samples <- names(specmine_2d_dataset$data)[which(names(specmine_2d_dataset$data)==spec_samples)]
    } else if (typeof(spec_samples)=="double") {
      samples <- names(specmine_2d_dataset$data)[spec_samples]
    }
  } else {
    if (!is.null(meta)) {
      snr_spectra <- snr_all(spect = specmine_2d_dataset, metadata = meta)
      samples <- rownames(snr_spectra[stats::ave(-snr_spectra$SNR,snr_spectra$Meta,FUN = rank) <= 2,])
    } else {
      snr_spectra <- snr_all(spect = specmine_2d_dataset)
      samples <- rownames(snr_spectra)[c(1:2,nrow(snr_spectra)-1,nrow(snr_spectra))]
    }
  }
  
  butons <- list()
  for (i in (1:(length(samples)+1))) {
    temp <- list()
    temp$method <- "restyle"
    if (i == length(samples)+1) {
      temp$args <- list("visible", c(rep(T, length.out = length(samples))))
      temp$label <- "All samples"
    } else if (i != length(samples)+1) {
      temp$args <- list("visible",c(rep(F,length.out = length(samples))))
      temp$args[[2]][i] <- T
      temp$label <- samples[i]}
    butons[[i]] <- temp
  }
  
  p1 <- plotly::plot_ly(showscale = F, hovertemplate = paste(   
    "F1: %{x:.0f}<br>",
    "F2: %{y:.0f}<br>",
    "Intensity: %{z:.0f}"))
  
  state <- specmine_2d_dataset$metadata[samples[1],meta]
  index <- 2
  for (sample in samples) {
    if (!is.null(meta)) {
      new_state <- specmine_2d_dataset$metadata[sample,meta]
      if (new_state != state) {
        state <- new_state
        index <- index + 1
      }
      colorstate <- clrs[index]
      p1 <- plotly::add_surface(p1, z = specmine_2d_dataset$data[[sample]], name = paste(state,sample,sep = ": "), visible = F, colorscale = colorstate, showlegend = T, legendgroup = state)
    } else {
      colorstate <- clrs[index]
      p1 <- plotly::add_surface(p1, z = specmine_2d_dataset$data[[sample]], name = sample, visible = F, showlegend = T, colorscale = colorstate, showlegend = T)
      index <- index + 1
    }
  }  
  
  
  initial_val <- as.numeric(rownames(specmine_2d_dataset$data[[samples[1]]]))[1]
  final_val <- as.numeric(rownames(specmine_2d_dataset$data[[samples[1]]]))[nrow(specmine_2d_dataset$data[[samples[1]]])]
  if (initial_val < 0){
    initial_val <- 0
  }
  if (final_val < 0){
    final_val <- 0
  }
  
  initial_val2 <- as.numeric(colnames(specmine_2d_dataset$data[[samples[1]]]))[1]
  if (initial_val2 < 0){
    initial_val2 <- 0
  }
  final_val2 <- as.numeric(colnames(specmine_2d_dataset$data[[samples[1]]]))[ncol(specmine_2d_dataset$data[[samples[1]]])]
  if (final_val2 < 0){
    final_val2 <- 0
  }
  
  seq1 <- seq(initial_val, final_val, length.out = 5)
  seq2 <- seq(initial_val2, final_val2, length.out = 5)
  
  seq_dim1 <- seq(0,nrow(specmine_2d_dataset$data[[samples[1]]]),length.out = 5)
  seq_dim2 <- seq(0,ncol(specmine_2d_dataset$data[[samples[1]]]),length.out = 5)
  
  xaxis_1 <- list(title = 'F2', tickvals = seq_dim2, ticktext = seq2)
  yaxis_2 <- list(title = 'F1', tickvals = seq_dim1, ticktext = seq1)
  
  p1 <- plotly::layout(p1, title = title_spectra,
                       scene = list(xaxis = xaxis_1, yaxis = yaxis_2, zaxis = list(title = 'Intensity')),
                       updatemenus = list(list(type = "dropdown", y = 0.8, buttons = butons)), legend = list(title=list(text=paste0('<b>',meta))))
  p1
}


##################################
#####PLOT OF MS AND NMR PEAKS#####
##################################

plot_peaks=function(dataset, column.class, samples = NULL, 
                    variable.bounds = NULL, xlab = NULL, ylab = NULL,
                    legend.place = "topright", cex = 0.8, reverse.x = FALSE, p.size=0.5, ...){
  
  if (is.null(xlab)) xlab = get_x_label(dataset)
  if (is.null(ylab)) ylab = get_value_label(dataset)
  
  if (is.null(variable.bounds)){
    variables = rownames(dataset$data)
    vars = variables
  } 
  else {
    x.vars = get_x_values_as_num(dataset)
    variables = rownames(dataset$data)[x.vars > variable.bounds[1] & x.vars < variable.bounds[2]] 
    vars = variables
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
  
  matplot(vars, dataset$data[variables,samples], col=as.integer(metadata),
          xlab=xlab, ylab=ylab, xlim=xlim, type="p", pch=16, cex=p.size, ...)
  
  if (legend.place != "none")
    legend(legend.place, levels(metadata), cex=cex, fill = sort(as.integer(factor(levels(metadata))))) 
}




##########################################################################################################
## Multiplot from ggplot2 - function taken from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_%28ggplot2%29/

multiplot <- function(plots, plotlist=NULL, file, cols=1, layout=NULL) {
  # Make a list from the ... arguments and plotlist
  plots <- c(plots, plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = grid::viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
