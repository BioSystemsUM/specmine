names = c("genefilter", "xcms", "impute")

for (n in names){
	if (!require(n)){
		source("http://bioconductor.org/biocLite.R")
		biocLite(n)
	}
}
