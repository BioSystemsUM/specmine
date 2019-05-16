globalVariables("group")

names = c("genefilter", "xcms", "impute", "MAIT")

for (n in 1:length(names)){
	if (!(names[n] %in% rownames(installed.packages()))){
		source("http://bioconductor.org/biocLite.R")
		biocLite(names[n])
	}
}

if (!"rcytoscapejs" %in% rownames(installed.packages())){
	if (!"devtools" %in% rownames(installed.packages())) {
	  if(!"memoise" %in% rownames(installed.packages())){install.packages("memoise", repos="https://cloud.r-project.org")}
	  install.packages("devtools", repos="https://cloud.r-project.org")
	}
	devtools::install_github('cytoscape/r-cytoscape.js@v0.0.7')
}

#if (!("mzmatch.R" %in% rownames(installed.packages()))){
#	source ("http://puma.ibls.gla.ac.uk/mzmatch.R/install_mzmatch.R")
#}
