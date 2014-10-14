names = c("genefilter", "xcms", "impute")

for (n in names){
	if (!require(n, character.only = TRUE)){
		source("http://bioconductor.org/biocLite.R")
		biocLite(n, character.only = TRUE)
	}
}
