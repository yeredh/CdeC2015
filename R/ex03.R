setwd("C:/Users/Yered/Dropbox/CdeC")

library(GEOquery)
gse <- getGEO("GSE21653", GSEMatrix=TRUE)
show(gse)

filePaths = getGEOSuppFiles('GSE21653')
filePaths
