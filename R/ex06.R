library("ALL")
data("ALL")
bcell = grep("^B", as.character(ALL$BT))
moltyp = which(as.character(ALL$mol.biol) %in% c("NEG", "BCR/ABL"))

ALL_bcrneg = ALL[, intersect(bcell, moltyp)]
ALL_bcrneg$mol.biol = factor(ALL_bcrneg$mol.biol)

library("hgu95av2.db")
library("genefilter")
ALLfilt_bcrneg = nsFilter(ALL_bcrneg, var.cutoff=0.5)$eset

library("GSEABase")

gsc = GeneSetCollection(ALLfilt_bcrneg,
                        setType=KEGGCollection())
Am = incidence(gsc)
dim(Am)

nsF = ALLfilt_bcrneg[colnames(Am),]


selectedRows = (rowSums(Am)>10)
Am2 = Am[selectedRows, ]

library("Category")
set.seed(123)
NPERM = 1000
pvals = gseattperm(nsF, nsF$mol.biol, Am2, NPERM)
pvalCut = 0.025
lowC = names(which(pvals[, 1]<=pvalCut))
highC = names(which(pvals[, 2]<=pvalCut))

library("KEGG.db")
head(getPathNames(lowC))
head(getPathNames(highC))
