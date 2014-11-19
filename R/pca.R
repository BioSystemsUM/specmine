############################################################################
################################PCA#########################################
############################################################################

# perform pca analysis - classical

pca.analysis.dataset = function(dataset, scale = T, center = T, 
                                write.file = F, file.out = "pca", ...) {
  
	pca.result = prcomp(t(dataset$data), center = center, scale. = scale, ...)
  if (write.file) {
    write.csv(pca.result$x, file=paste(file.out,"_scores.csv",sep=""))
    write.csv(pca.result$rotation, file=paste(file.out,"_loadings.csv", sep= ""))
  }
	pca.result
}

# returns information about importance of the PC's
# pcs - PCs to get; sd - get std dev; prop - get proportion of variance; cumul - get cumulative
# min.cum - allows to define minimum cumulative % of variance
pca.importance = function(pca.res, pcs = 1:length(pca.res$sdev), sd = T, prop = T, cumul = T, min.cum = NULL)
{
  rows = c()
  if (sd) rows = c(1)
  if (prop) rows = c(rows, 2)
  if (cumul) rows = c(rows, 3)
  
  if (class(pca.res) == "prcomp"){
	  if (!is.null(min.cum)) {
		cum = summary(pca.res)$importance[3,]
		pcs = 1:(min(which(cum > min.cum)))
	  }

	  res = summary(pca.res)$importance
  } else if (class(pca.res) == "princomp"){
	  vars = pca.res$sdev^2
	  vars = vars/sum(vars)
	  cum = cumsum(vars)
	  if (!is.null(min.cum)) {
		pcs = 1:(min(which(cum > min.cum)))
	  }
	  res = rbind("Standard deviation" = pca.res$sdev, "Proportion of Variance" = vars, "Cumulative Proportion" = cum)
  }
  res[rows, pcs]
}

# robust PCA analysis 

# center - how the data will be centered "mean" or "median" (or NULL if nore)
# scale - how the data will be scaled "sd" or "mad" (or NULL if none)
# k - number of PCs to compute

# returns objects of class princomp
pca.robust = function(dataset, center = "median", scale = "mad", k = 10,
                      write.file = F, file.out = "robpca", ...)
{
  require(pcaPP)
  pca.res = PCAgrid(t(dataset$data), k = k, center = center, scale = scale, scores = T, ...)
  if (write.file) {
    write.csv(pca.result$scores, file=paste(file.out,"_scores.csv",sep=""))
    write.csv(pca.result$loadings, file=paste(file.out,"_loadings.csv", sep= ""))
  }
  pca.res
}


########################## PCA PLOTS ##################################

#scree plot
pca.screeplot = function(pca.result, num.pcs = NULL, cex.leg = 0.8, leg.pos = "right", 
                         lab.text = c("individual percent","cumulative percent"), 
                         fill.col = c("blue","red"), ylab = "Percentage", xlab = "Principal components",
                         ...){
  importance = pca.importance(pca.result)
  if (is.null(num.pcs)) num.pcs = dim(importance)[2]
  par(mfrow=c(1,1))
  matplot(seq(1, num.pcs), data.frame(t(importance[2:3,1:num.pcs]*100)), type="l", lty=1, col=fill.col, 
          xaxt='n', ylab= ylab,xlab=xlab, ...)
  legend(leg.pos, lab.text, cex=cex.leg, fill=fill.col)
  axis(1, at = 1:dim(importance)[2],labels= colnames(importance))
}

#2d scores plot
pca.scoresplot2D = function(dataset, pca.result, column.class, pcas = c(1,2), labels = FALSE, 
                            ellipses = FALSE, pallette = 2)
{
  require(ggplot2)
  if (class(pca.result) == "prcomp"){
	scores = pca.result$x
  } else if (class(pca.result) == "princomp"){
	scores = pca.result$scores
  }
  pca.points = data.frame(scores[,pcas])
  names(pca.points) = c("x","y")
  pca.points$group = dataset$metadata[,column.class]
  pca.points$label = colnames(dataset$data)
  pca.plot = ggplot(data = pca.points, aes(x=x, y=y,colour=group)) + geom_point(size=3, alpha=1) +
    scale_colour_brewer(type = "qual", palette=pallette) + xlab(paste("PC",pcas[1],sep="")) + ylab(paste("PC",pcas[2],sep=""))
  if (labels){
    pca.plot = pca.plot + geom_text(data = pca.points, aes(x,y,label=label),hjust=-0.1, vjust=0)
  }
  if (ellipses){
    df.ellipses = calculate.ellipses(pca.points)
    pca.plot = pca.plot + geom_path(data=df.ellipses, aes(x=x, y=y,colour=group), size=1, linetype=2) 
  }
  pca.plot
}

#3d scores plot
pca.scoresplot3D.rgl = function(dataset, pca.result, column.class, pcas = c(1,2,3), size = 1, 
                            labels = FALSE) {
  require(rgl)
  if (class(pca.result) == "prcomp"){
	scores = pca.result$x
  } else if (class(pca.result) == "princomp"){
	scores = pca.result$scores
  }
  plot3d(scores[,pcas], type = "s", col = as.integer(dataset$metadata[,column.class]),
         size=size)
  if (labels){
    text3d(scores[,pcas],texts=colnames(dataset$data), cex=0.6)
  }
}

pca.scoresplot3D = function(dataset, pca.result, column.class, pcas=c(1,2,3))
{
  require(scatterplot3d)
  if (class(pca.result) == "prcomp"){
	scores = pca.result$x
  } else if (class(pca.result) == "princomp"){
	scores = pca.result$scores
  }
  classes = dataset$metadata[,column.class]
  scatterplot3d(scores[,pcas], color=as.integer(dataset$metadata[,column.class]), pch=17)
  legend(-1.5, 2.5, levels(classes), col = 1:length(classes), cex = 0.7, pt.cex = 1, pch= 17)
}

#biplots
pca.biplot= function(pca.result, cex = 0.8, ...) {
  biplot(pca.result, cex = cex, ...)
}

pca.biplot3D = function(dataset, pca.result, column.class, pcas = c(1,2,3)){
  if (class(pca.result) == "prcomp"){
	scores = pca.result$x
	rotation = pca.result$rotation
  } else if (class(pca.result) == "princomp"){
	scores = pca.result$scores
	rotation = pca.result$loadings
  }  
  pca.scoresplot3D.rgl(dataset, pca.result, column.class, pcas)
  text3d(scores[,pcas], texts=colnames(dataset$data), cex=0.6)
  text3d(rotation[,pcas], texts = rownames(rotation), col = "red", cex=0.6)
  coords = NULL
  for (i in 1:nrow(rotation)){
    coords = rbind(coords, rbind(c(0,0,0), rotation[i, pcas]))
  }
  lines3d(coords, col="red", lwd = 4)
}

#pca pairs plot
pca.pairs.plot = function(dataset, column.class, pca.result, pcas = c(1,2,3,4,5), ...){
  require(GGally)
  if (class(pca.result) == "prcomp"){
	scores = pca.result$x
  } else if (class(pca.result) == "princomp"){
	scores = pca.result$scores
  }  
  pairs.df = data.frame(scores[,pcas])
  pairs.df$group = dataset$metadata[,column.class]
  ggpairs(pairs.df, colour = 'group', ...)
}

#kmeans clustering with 3 PCs
pca.kmeans.plot3D = function(dataset, pca.result, num.clusters = 3, pcas = c(1,2,3), 
                             kmeans.result = NULL, labels = FALSE, size = 1,...) {
  require(rgl)
  if (is.null(kmeans.result)){
    kmeans.result = clustering(dataset, method = "kmeans", num.clusters = num.clusters)
  }
  plot3d(pca.result$x[,pcas], type = "s", col = kmeans.result$cluster, size=size,...)
  if (labels){
    text3d(pca.result$x[,pcas],adj = c(1.2,1.2), texts=colnames(dataset$data), cex=0.6)
  }
}


#kmeans clustering with 2 first PCs
pca.kmeans.plot2D = function(dataset, pca.result, num.clusters = 3, pcas = c(1,2), 
                             kmeans.result = NULL, labels = FALSE, ellipses = FALSE){
  require(ggplot2)
  if (is.null(kmeans.result)){
    kmeans.result = clustering(dataset, method = "kmeans", num.clusters = num.clusters)
  }
  pca.points = data.frame(pca.result$x[,pcas])
  names(pca.points) = c("x","y")
  pca.points$group = factor(kmeans.result$cluster)
  pca.points$label = colnames(dataset$data)
  pca.plot = ggplot(data = pca.points, aes(x=x, y=y,colour=group)) + geom_point(size=3, alpha=.6) +
    scale_colour_brewer(palette="Set1") + xlab(paste("PC",pcas[1],sep="")) + ylab(paste("PC",pcas[2],sep=""))
  if (labels){
    pca.plot = pca.plot + geom_text(data = pca.points, aes(x,y,label=label),hjust=-0.1, vjust=0, size = 3)
  }
  if (ellipses){
    df.ellipses = calculate.ellipses(pca.points)
    pca.plot = pca.plot + geom_path(data=df.ellipses, aes(x=x, y=y,colour=group), size=1, linetype=2) 
  }
  pca.plot
}


#pca pairs with kmeans clusters plot
pca.pairs.kmeans.plot = function(dataset, pca.result, num.clusters = 3, kmeans.result = NULL, pcas = c(1,2,3,4,5)){
  require(GGally)
  if (is.null(kmeans.result)){
    kmeans.result = clustering(dataset, method = "kmeans", num.clusters = num.clusters)
  }
  pairs.df = data.frame(pca.result$x[,pcas])
  pairs.df$group = factor(kmeans.result$cluster)
  ggpairs(pairs.df, colour = 'group')
}

#draw ellipses
calculate.ellipses = function(data){
  require(ellipse)
  df_ell <- data.frame()
  for(g in levels(data$group)){
    df_ell <- rbind(df_ell, cbind(as.data.frame(with(data[data$group==g,], ellipse(cor(x, y), 
                      scale=c(sd(x),sd(y)), centre=c(mean(x),mean(y))))),group=g))
  }
  df_ell
}



