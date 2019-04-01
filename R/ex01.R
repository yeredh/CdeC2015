# dagdata
library(dagdata)
data(maPooling)

setwd("C:/Users/Yered/Dropbox/CdeC")

library(Biobase)


# RData
load("maPooling.RData")
e <- maPooling
head(rownames(e))


# Plataforma
annotation(e)

# Paquetes para anotacion
library(rae230a.db)
library(AnnotationDbi)

# Keys: campos que se pueden como palabras claves 
keytypes(rae230a.db)
# Por ejemplos, podemos usar los nombres de las sondas para acceder a los otros campos de la anotacion
head(keys(rae230a.db, keytype="PROBEID"))
# En nuestas muestras
head(rownames(e))

# Nombres de los Genes
# Ensembl (European Bioinformatics Institute \& Wellcome Trust Sanger Institute )
# Entrez ( National Center for Biotechnology Information)
# Symbol ( Human Genome Organisation)
res <- select(rae230a.db, keys=rownames(e),
              columns=c("ENTREZID","ENSEMBL","SYMBOL"), 
              keytype="PROBEID")

head(res)

# Vamos a agregar esta informacion a nuestras muestras
idx <- match(rownames(e), res$PROBEID)
head(rownames(e))
head(res$PROBEID,7)

fData(e) <- res[idx,]
head(fData(e),10)

# Vamos a asegurarnos que los nombres corresponden
all.equal(fData(e)$PROBEID, rownames(e))
