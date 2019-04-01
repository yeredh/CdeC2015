# leer datos procesados
eset = readRDS("Alz.RDS")

# grupos
table(eset$Alz)

# Seleccionar grupos
types = c("Control", "Severe")
ind = which(as.character(eset$Alz) %in% types)
Alz_severe = eset[,ind]
Alz_severe$Alz = factor(Alz_severe$Alz)
# eliminar valores faltantes
Alz_severe <- Alz_severe[rowSums(is.na(exprs(Alz_severe)))==0,]

# filtrar sondas por variabilidad
library("genefilter")
sds = rowSds(exprs(Alz_severe))
# eliminar sondas con nula variabilidad
Alz_severe <- Alz_severe[!is.na(sds),]
sds = sds[!is.na(sds)]
sh = shorth(sds,na.rm=TRUE)
sh

# histograma de las desviaciones estandar de las sondas
hist(sds, breaks=50, col="mistyrose",
     xlab="desviación estándar",main="")
abline(v=sh, col="blue", lwd=3, lty=2)

# Seleccionar sondas con desveiacion estandar por arriba de sh
Alz_filt = Alz_severe[sds>=sh, ]
dim(exprs(Alz_filt))
dim(exprs(Alz_severe))


table(Alz_severe$Alz)
# Prueba de t
tt = rowttests(Alz_filt, "Alz")
names(tt)
# Histograma de los valores p
hist(tt$p.value, breaks=50, col="mistyrose", xlab="valor p",main="")


# Ajustar valores p
library("multtest")
mt = mt.rawp2adjp(tt$p.value, proc="BH")
# Obtener nombre de las 10 sondas mas significativas
g = featureNames(Alz_filt)[mt$index[1:10]]


# Obtener nombres de los genes
library("hgu133plus2.db")
res_symbols =links(hgu133plus2SYMBOL[g])


# terminos funcionales
library(GO.db)
goterms = unlist(Term(GOTERM))
go_res = links(hgu133plus2GO[g])

res <- c()
i = 1
while(i <= length(go_res$go_id)){
  tmp = cbind(go_res$probe_id[i],goterms[names(goterms) %in% go_res$go_id[i]])
  res <- rbind(res,tmp)
  i=i+1
}
rownames(res)= NULL
colnames(res)=c("probe.id","GO.term")
terminos_funcionales = merge(res_symbols,res,by.y="probe.id",by.x="probe_id")
terminos_funcionales


# GSEA
library("GSEABase")
# Usar grupos de genes de KEGG
gsc = GeneSetCollection(Alz_filt,setType=KEGGCollection())
Am = incidence(gsc)
dim(Am)
# usar solo grupos de genes con mas de 10 elementos
nsF = Alz_filt[colnames(Am),]
selectedRows = (rowSums(Am)>10)
Am2 = Am[selectedRows, ]
# obtener valores p de GSEA
library("Category")
set.seed(123)
NPERM = 1000
pvals = gseattperm(nsF, nsF$Alz, Am2, NPERM)
pvalCut = 0.025
lowC = names(which(pvals[, 1]<=pvalCut))
highC = names(which(pvals[, 2]<=pvalCut))
# obtener nombres de los grupos de genes
library("KEGG.db")
head(getPathNames(lowC))
head(getPathNames(highC))