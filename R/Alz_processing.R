

library(Biobase)
library(GEOquery)

# Leer archivo SOFT descargado de GEO
gse <- getGEO(filename="GSE28146_family.soft.gz")

# Nombres de las sondas
probesets <- Table(GPLList(gse)[[1]])$ID
# Obtener los valores de para cada muestra
data.matrix <- do.call('cbind',lapply(GSMList(gse),function(x){tab <- Table(x)
                                                               mymatch <- match(probesets,tab$ID_REF)
                                                               return(tab$VALUE[mymatch])}))
data.matrix <- apply(data.matrix,2,function(x) {as.numeric(as.character(x))})
data.matrix <- log2(data.matrix)
data.matrix[1:5,]

# Crear objecto de ExpressionSet
rownames(data.matrix) <- probesets
colnames(data.matrix) <- names(GSMList(gse))
pdata <- data.frame(samples=names(GSMList(gse)))
rownames(pdata) <- names(GSMList(gse))
pheno <- as(pdata,"AnnotatedDataFrame")
eset <- new('ExpressionSet',exprs=data.matrix,phenoData=pheno)
eset

# Agregar dato de anotation
# plataforma de microarreglo
annotation(eset) <- "hgu133plus2"

# Status de Enfermedad
tmp=t(sapply(GSMList(gse),function(x) {Meta(x)$characteristics_ch1}))
eset$Alz <- gsub("disease status: ","",tmp[,3])

# guardar datos procesados
saveRDS(eset,"Alz.RDS")

