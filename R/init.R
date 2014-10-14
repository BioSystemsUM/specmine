names = c("genefilter", "xcms", "impute")

for (n in 1:length(names)){
	if (!(names[n] %in% rownames(installed.packages()))){
		source("http://bioconductor.org/biocLite.R")
		biocLite(names[n])
	}
}
