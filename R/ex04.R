library("Biobase")
library("ALL")
library("genefilter")
data("ALL")

ALL$BT 
bcell = grep("^B", as.character(ALL$BT))

types = c("NEG", "BCR/ABL")
moltyp = which(as.character(ALL$mol.biol) %in% types)

ALL_bcrneg = ALL[, intersect(bcell, moltyp)]
ALL_bcrneg$mol.biol = factor(ALL_bcrneg$mol.biol)


library("genefilter")
sds = rowSds(exprs(ALL_bcrneg))
sh = shorth(sds)
sh

hist(sds, breaks=50, col="mistyrose",
     xlab="desviación estándar",main="")
abline(v=sh, col="blue", lwd=3, lty=2)


ALLsfilt = ALL_bcrneg[sds>=sh, ]
dim(exprs(ALLsfilt))
dim(exprs(ALL_bcrneg))


table(ALLsfilt$mol.biol)
tt = rowttests(ALLsfilt, "mol.biol")
names(tt)

hist(tt$p.value, breaks=50, col="mistyrose", xlab="valor p",main="")


library("multtest")
mt = mt.rawp2adjp(tt$p.value, proc="BH")


g = featureNames(ALLsfilt)[mt$index[1:10]]


library("hgu95av2.db")
links(hgu95av2SYMBOL[g])


mb = ALLsfilt$mol.biol
y = exprs(ALLsfilt)[g[1],]
ord = order(mb)
plot(y[ord], pch=c(1,16)[mb[ord]],
     col=c("black", "red")[mb[ord]],
     main=hgu95av2SYMBOL[g[1]], ylab="Expresión",
     xlab="Muestras")

