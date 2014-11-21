########################################################################
######################## LOW LEVEL FUSION ##############################
########################################################################

"low.level.fusion" = function(datasets){
	sample.names = colnames(datasets[[1]]$data)
	for (i in 2:length(datasets)){
		sample.names = intersect(sample.names, colnames(datasets[[i]]$data))
	}
	subsets = list()
	for (i in 1:length(datasets)){
		if (!is.null(datasets[[i]]$metadata)) r.factors = T
		else r.factors = F
		subsets[[i]] = subset.samples(datasets[[i]], samples = sample.names, rebuild.factors = r.factors)
	}
	
	ds.fused = fusion.merge(subsets)
	ds.fused
}

# assuming rownames are not duplicated and first dataset contains all the metadata 
"fusion.merge" = function(datasets){
	ds.fused = datasets[[1]]
	for (i in 2:length(datasets)){
		ds.fused$data = rbind(ds.fused$data, datasets[[i]]$data)
	}
	ds.fused$metadata = datasets[[1]]$metadata
	ds.fused
}
