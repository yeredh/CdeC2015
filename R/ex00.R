setwd("C:/Users/Yered/Dropbox/CdeC")

# Working Directory
wd <- getwd()
basedir <- paste0(wd, "/celfiles")
setwd(basedir)

# Load sample information
library(affy)
tab <- read.delim("sampleinfo.txt",check.names=FALSE,as.is=TRUE)
rownames(tab) <- tab$filenames
tab

# List raw files
fns <- list.celfiles()
fns
fns %in% tab[,1]

# read in CEL files
ab <- ReadAffy(phenoData=tab)
dim(pData(ab))

# Platform
annotation(ab)

# Processing
e <- rma(ab)


# Dettach affy package to load oligo package
detach("package:affy")
library(oligo)

# set wd and load sample info
basedir <- paste0(wd,"/celfiles")
setwd(basedir)
tab <- read.delim("sampleinfo.txt",check.names=FALSE,as.is=TRUE)

# list files
fns <- list.celfiles(listGzipped=TRUE)
fns %in% tab[,1]

# read & process data
pd <- as(tab, "AnnotatedDataFrame")
efs <- read.celfiles(filenames=tab[,1],phenoData=pd,sampleNames=sampleNames(pd))
setwd(wd)
e <- rma(efs)

