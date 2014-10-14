names = c("genefilter", "xcms", "impute")

for (n in names){
	if (!(n %in% rownames(installed.packages()))){
		source("http://bioconductor.org/biocLite.R")
		biocLite(n, character.only = TRUE)
	}
}
