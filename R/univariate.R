# univariate tests
# type: "anova", "ttest", "foldchange" or "correlation" analysis
# column.class: which column of metadata will be used (index or name)
# ref.class: reference value of metadata variable that will be used in fold change
# method: method used in correlation analysis

### TO REMOVE THIS HIGH-LEVEL FUNCTION - DOESN'T MAKE SENSE ###
"univariate.analysis" = function(dataset, type = "anova", column.class, ref.class, 
                                 method = "pearson"){
	if (type == "anova"){
		result = aov.all.vars(dataset, column.class)
	} 
  else if (type == "ttest"){
		result = tTests.dataset(dataset, column.class)
	} 
  else if (type == "foldchange"){
		result = fold.change(dataset, column.class, ref.class)
	} 
  else if (type == "correlation"){
		result = correlations.dataset(dataset, method)
	}
  else stop("Type of univariate analysis not defined")
	result
}


# ANALYSIS OF VARIANCE

"aov.one.var" = function(dataset, x.val, groups, doTukey= T)
{
  values = get.data.values(dataset, x.val)
  resaov = aov(values ~ groups)
  res = list()
  res$pvalue = summary(resaov)[[1]]$'Pr(>F)'[1]
  res$log = -1*log10(res$pvalue)
  res$fdr = p.adjust(res$pvalue, "fdr")
  if (doTukey) {
    resTukey = TukeyHSD(resaov)
    inx <- resTukey$groups[,"p adj"] <= 0.05
    res$tukey = paste(rownames(resTukey$groups)[inx], collapse="; ")
  }
  res
}

"aov.all.vars" = function(dataset, column.class, doTukey= T, 
                          write.file = F, file.out = "anova-res.csv" )
{
  groups = dataset$metadata[,column.class]
  pvalues = c()
  logs = c()
  fdr = c()
  if (doTukey) tukeys = c()
  for(i in 1:nrow(dataset$data))
  {
    resi = aov.one.var(dataset, rownames(dataset$data)[i], groups, doTukey)
    pvalues[i] = resi$pvalue
    logs[i] = resi$log
    if (doTukey) tukeys[i] = resi$tukey
  }
  fdr = p.adjust(pvalues,"fdr")
  res = data.frame(pvalues, logs, fdr)
  if (doTukey) res$tukey = tukeys
  rownames(res) = rownames(dataset$data)
  
  aov.table = res[order(pvalues),]
  if (write.file) write.csv(aov.table, file=file.out)
  aov.table
}

"multifactor.aov.onevar" = function(dataset, x.val, metadata.vars, combination) {
  values = get.data.values(dataset, x.val)
  sub.df = dataset$metadata[,metadata.vars]
  sub.df = cbind(values,sub.df)
  terms = names(sub.df)[2:ncol(sub.df)]
  terms = cbind(terms, combination)
  resaov = aov(reformulate(terms, "values"), data = sub.df)
  aov.summary = summary(resaov)[[1]]
  aov.summary
}

"multifactor.aov.all.vars" = function(dataset, metadata.vars, combination)
{
  #m = matrix(NA, nrow(dataset$data), length(metadata.vars))
  #rownames(m) = rownames(dataset$data)
  m = vector("list",nrow(dataset$data))

  for(i in 1:nrow(dataset$data))
  {
    m[[i]] = multifactor.aov.onevar(dataset, rownames(dataset$data)[i], metadata.vars, combination)
  }
  names(m) = rownames(dataset$data)

  #aov.table = as.data.frame(m)
  #if (is.numeric(metadata.vars)) colnames(aov.table) = colnames(dataset$metadata)[metadata.vars]
  #else colnames(aov.table) = metadata.vars
  #if (write.file) write.csv(aov.table, file=file.out)
  #aov.table
  m
}

multifactor.aov.pvalues.table = function(multifactor.aov.results, write.file = F, file.out = "multi-anova-pvalues.csv"){
	num_vars = length(multifactor.aov.results[[1]]$'Pr(>F)') - 1
	m = matrix(NA, length(multifactor.aov.results), num_vars)
	rownames(m) = names(multifactor.aov.results)
	for (i in 1:length(multifactor.aov.results)){
		m[i,] = multifactor.aov.results[[i]]$'Pr(>F'[1:num_vars]
	}
	aov.table = as.data.frame(m)
	colnames(aov.table) = trim(head(rownames(multifactor.aov.results[[1]]),-1))
	if (write.file) write.csv(aov.table, file = file.out)
	aov.table	
}

multifactor.aov.varexp.table = function(multifactor.aov.results, write.file = F, file.out = "multi-anova-varexp.csv"){
  num_vars = length(multifactor.aov.results[[1]]$'Sum Sq') # - 1
  m = matrix(NA, length(multifactor.aov.results), num_vars)
  rownames(m) = names(multifactor.aov.results)
  for (i in 1:length(multifactor.aov.results)){
    ssq = sum(multifactor.aov.results[[i]]$'Sum Sq')
    m[i,] = multifactor.aov.results[[i]]$'Sum Sq'[1:num_vars] /ssq
  }
  aov.table = as.data.frame(m)
  colnames(aov.table) = trim(rownames(multifactor.aov.results[[1]]))
  if (write.file) write.csv(aov.table, file = file.out)
  aov.table	
}
trim <- function( x ) {
  gsub("(^[[:space:]]+|[[:space:]]+$)", "", x)
}

plot.anova = function(dataset, anova.results, anova.threshold = 0.01, reverse.x = F) {
  orig.ord = intersect (get.x.values.as.text(dataset), rownames(anova.results))
  anova.orig = anova.results[orig.ord,]
  anova.lower = which(anova.orig$pvalues < anova.threshold)

  cols = vector("character", nrow(anova.orig))
  for(i in 1:nrow(anova.orig))
    if (i %in% anova.lower) cols[i] = "blue"
    else cols[i] = "gray"
  
  vars = rownames(anova.orig)
  if (reverse.x){
	xlim = c(max(as.numeric(vars)), min(as.numeric(vars)))
  } else {
	xlim = range(as.numeric(vars))
  }

  plot(vars,anova.orig$"logs", xlab = get.x.label(dataset), ylab = "-log10(p)", col = cols, pch = 19, xlim = xlim)
  #axis(1, at = 1:length(vars),labels = vars, xlim = xlim)
  abline(h = -log10(anova.threshold), col = "lightblue")
}

##################### FOLD CHANGE ############################
fold.change.var = function(dataset, metadata.var, variables, threshold.min.fc = NULL, 
								write.file = F, file.out = "fold_change_reverse.csv"){
	
	samp.classes = dataset$metadata[,metadata.var]
	datamat = t(dataset$data)
	means = matrix(ncol = 2, nrow = length(levels(samp.classes)))
	for (i in 1:length(levels(samp.classes))){
		means[i,1] = mean(datamat[which(samp.classes == levels(samp.classes)[i]), variables[1]])
		means[i,2] = mean(datamat[which(samp.classes == levels(samp.classes)[i]), variables[2]])
	}	
	#indexes = as.integer(samp.classes)
	#dataset = aggregate.samples(dataset, indexes, aggreg.fn = "mean")
	#datamat = t(dataset$data)
	#mean1 = rowMeans(datamat)
	
	fc.all = means[,1] / means[,2]
	fc.log = log2(fc.all)
	fc.res = data.frame(fc.all, fc.log)
	rownames(fc.res) = levels(samp.classes)
	colnames(fc.res) = c("FoldChange", "log2(FC)")
	  if (!is.null(threshold.min.fc)) {
		fc.res = subset(fc.res, FoldChange > threshold.min.fc | FoldChange < 1/threshold.min.fc)
	  }
	  fold.order = order(abs(fc.res[,2]), decreasing = T)
	  fc.res = fc.res[fold.order,,drop=F]
	if (write.file) write.csv(fc.res, file = file.out)
	fc.res 
}

fold.change = function(dataset, metadata.var, ref.value, threshold.min.fc = NULL,
                       write.file = F, file.out = "fold_change.csv" ) {
  datamat = dataset$data
  samp.classes = dataset$metadata[,metadata.var]
	mean1 = rowMeans(datamat[,which(samp.classes == ref.value)])
	mean2 = rowMeans(datamat[,which(samp.classes != ref.value)])
	fc.all = mean2/mean1
	fc.log = log2(fc.all)
	fc.res = data.frame(fc.all, fc.log)
	rownames(fc.res) = rownames(datamat)
	colnames(fc.res) = c("FoldChange", "log2(FC)")
  if (!is.null(threshold.min.fc)) {
    fc.res = subset(fc.res, FoldChange > threshold.min.fc | FoldChange < 1/threshold.min.fc)
  }
  fold.order = order(abs(fc.res[,2]), decreasing = T)
  fc.res = fc.res[fold.order,,drop=F]
  if (write.file) write.csv(fc.res, file = file.out)
	fc.res 
}

plot.fold.change = function(dataset, fc.results, fc.threshold, plot.log = T, var = F, xlab = "") {
  if (var == F){
	orig.ord = intersect (get.x.values.as.text(dataset), rownames(fc.results))
	fc.orig = fc.results[orig.ord,]
	xlabel = get.x.label(dataset)
  } else {
	fc.orig = fc.results
	xlabel = xlab
  }
  fc.higher = which(fc.orig$FoldChange > fc.threshold | 
                    fc.orig$FoldChange < 1/fc.threshold)
  cols = vector("character", nrow(fc.orig))
  for(i in 1:nrow(fc.orig))
    if (i %in% fc.higher) cols[i] = "blue"
    else cols[i] = "gray"
  if (plot.log) {
    max = max ( max(abs(fc.orig$"log2(FC)")), abs(log2(fc.threshold)) )
    plot(fc.orig$"log2(FC)", xlab = xlabel, ylab = "Log2FoldChange", 
         col = cols, pch = 19, ylim = c(-max,max), xaxt="n")
    axis(1, at = 1:length(rownames(fc.orig)),labels = rownames(fc.orig))
    abline(h = log2(fc.threshold), col = "lightblue")
    abline(h = -log2(fc.threshold), col = "lightblue")
    abline(h = 0)
  }
  else {
    max = max ( max(fc.orig$FoldChange), fc.threshold )
    min = min ( min(fc.orig$FoldChange), 1/fc.threshold) 
    plot(fc.orig$FoldChange, xlab = xlabel, ylab = "FoldChange", 
         col = cols, pch = 19, ylim = c(min,max), xaxt="n")
    axis(1, at = 1:length(rownames(fc.orig)),labels = rownames(fc.orig))
    abline(h = fc.threshold, col = "lightblue")
    abline(h = 1/fc.threshold, col = "lightblue")
    abline(h = 1)
  }
  
}

####################### T-TESTS ##################################

tTests.pvalue = function(datamat, samp.classes) {
	require(genefilter)
	p.value = try(rowttests(datamat, samp.classes)$p.value);
	if (class(p.value) == "try-error"){
		p.value = NA
	} else {
		names(p.value) = rownames(datamat)
	}
	res = data.frame(p.value)
  rownames(res) = names(p.value)
  res[order(p.value),,drop=F]
}

tTests.dataset = function(dataset, metadata.var, threshold = NULL,
                          write.file= F, file.out = "ttests.csv") {
	datamat = dataset$data
  samp.classes = dataset$metadata[,metadata.var]
  res.p.values = tTests.pvalue (datamat, samp.classes)
	p.log = -log10(res.p.values$p.value)
  fdr.p = p.adjust(res.p.values$p.value, "fdr")
	ttests.table = cbind(res.p.values, p.log, fdr.p)
	colnames(ttests.table) = c("p.value","-log10","fdr")
	if (!is.null(threshold)) {
    ttests.table = subset(ttests.table, p.value <= threshold)    
	}
	ttests.order = order(ttests.table[,1])
	ttests.table = ttests.table[ttests.order,,drop=F]
	if (write.file) write.csv(ttests.table, file=file.out)
	ttests.table
}

plot.ttests = function(dataset, tt.results, tt.threshold = 0.01) {
  orig.ord = intersect (get.x.values.as.text(dataset), rownames(tt.results))
  tt.orig = tt.results[orig.ord,]
  tt.lower = which(tt.orig$p.value < tt.threshold)

  cols = vector("character", nrow(tt.orig))
  for(i in 1:nrow(tt.orig))
    if (i %in% tt.lower) cols[i] = "blue"
    else cols[i] = "gray"
  
  plot(tt.orig$"-log10", xlab = get.x.label(dataset), ylab = "-log10(p)", col = cols, pch = 19, xaxt="n")
  axis(1, at = 1:length(rownames(tt.orig)),labels = rownames(tt.orig))
  abline(h = -log10(tt.threshold), col = "lightblue")
}

####################### VOLCANO PLOT FOR T-TESTS & FOLD CHANGES #######

volcano.plot.fc.tt = function(dataset, fc.results, tt.results, 
                              fc.threshold = 2, tt.threshold = 0.01) 
{
  ## test if compatible
  if (nrow(fc.results) != nrow(tt.results))
    stop("Fold change and ttest results are not compatible")
  int = intersect(rownames(fc.results), rownames(tt.results))
  if (length(int) != nrow(fc.results))
    stop("Fold change and ttest results are not compatible")
 
  orig.ord = intersect (get.x.values.as.text(dataset), rownames(tt.results))
  tt.orig = tt.results[orig.ord,]
  fc.orig = fc.results[orig.ord,]
  
  fc.new.values = sapply(fc.orig$FoldChange, function(x){ if (x < 1) return (1/x) else return (x) })
  
  to.color = which(tt.orig$p.value < tt.threshold & 
                  fc.new.values > fc.threshold)
  
  cols = vector("character", nrow(tt.orig))
  texts = vector("character", nrow(tt.orig))
  for(i in 1:nrow(tt.orig))
    if (i %in% to.color) {
      cols[i] = "blue"
      texts[i] = orig.ord[i]
    }
  else {
    cols[i] = "gray"
    texts[i] = ""
  }
  
  plot(fc.orig$"log2(FC)", tt.orig$"-log10", xlab = "log2(FC)", ylab = "-log10(p)", 
       col = cols, pch = 19)
  text(fc.orig$"log2(FC)", tt.orig$"-log10", texts, cex = 0.6, col = "blue", srt = -30, pos = 1)
  abline(h = -log10(tt.threshold), col = "lightblue")
  abline(v = log2(fc.threshold), , col = "lightblue")
  abline(v = -(log2(fc.threshold)), , col = 'lightblue') 
  rownames(dataset$data[to.color,])
}


####################### CORRELATIONS ############################
		
# method: pearson, kendall or spearman
# by.var - if T, correlations of variables (rows); if F, correlations of samples (columns)
correlations.dataset = function(dataset, method = "pearson", by.var = T) {
  
  if (by.var) data.to.cor = t(dataset$data)
  else data.to.cor = dataset$data
  
	cor.matrix = cor(data.to.cor, method = method)
  cor.matrix
}

heatmap.correlations = function(correlations, col = NULL, ...) {
  require(RColorBrewer)
  if (is.null(col)) colors = rev(colorRampPalette(brewer.pal(10, "RdBu"))(256))
  else colors = col
  heatmap(correlations, col = colors, ...)
}


correlation.test = function(dataset, x,y, method = "pearson", alternative = "two.sided", by.var = T){
	if (by.var) data.to.cor = data.frame(t(dataset$data))
	else data.to.cor = data.frame(dataset$data)

	result.cor = with(data.to.cor, cor.test(eval(as.name(x)),eval(as.name(y)),method = method,alternative = alternative))

	result.cor
}  


correlations.test = function(dataset, method = "pearson", by.var = T, alternative = "two.sided") {
	if (by.var) data.to.cor = t(dataset$data)
	else data.to.cor = dataset$data
	
	data.names = colnames(data.to.cor)
	cor.matrix = matrix(nrow=length(data.names)^2, ncol = 4)
	i = 1
	for (name in data.names){
		for (name2 in data.names){
			result.cor = with(data.frame(data.to.cor), correlation.test(dataset, as.name(name), as.name(name2), method = method, alternative = alternative, by.var = by.var))
			cor.estimate = result.cor$estimate
			cor.pvalue = result.cor$p.value
			cor.matrix[i,] = c(name, name2, cor.estimate, cor.pvalue)
			i = i+1
		}
	}
	colnames(cor.matrix) = c("i","j","cor","p.value")
	cor.matrix
}
