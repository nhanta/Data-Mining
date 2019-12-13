# 1. Install frma package
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("frma")
BiocManager::install("hgu133plus2frmavecs")
BiocManager::install("biomaRt")
BiocManager::install("GEOquery")

  # Use help for frma
browseVignettes("frma")

# 2. Load library

library(frma)
library(affy)
library(oligo)
library(hgu133plus2frmavecs)
require("biomaRt")
require(GEOquery)
require(Biobase)


# 3. Read CEL file
mypath <-"D:/Data Science/Data Science Core/Documents/Data Mining/Project/Ovary/C/"
data = ReadAffy(celfile.path = mypath, compress = TRUE)

# 4. Load an example AffyBatch

object = frma(data)
e <- exprs(object)
GNUSE(object,type="stats")

# 5. z-score, z-scores (s) and barcode

z_score <- barcode(object, output="z-score")
s = pnorm(z_score)
barcode <- barcode(object, output="binary")

n_samples = dim(z_score)[2]
for (i in (1:n_samples)){
  D = data.frame(z_score= z_score[,i], s=s[,i], barcode=barcode[,i])
  write.csv(D, file = paste0(mypath, i, ".csv"), row.names = TRUE)
  }

# 6. Annotate Affymetrix probesets to Gene symbols

  # 6.1. Get ids of gene
ids <- geneNames(data)

  # 6.2. Get table of gene information include gene ids and gene names

mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)
annotLookup <- getBM(
  mart=mart,
  attributes=c(
     "affy_hg_u133_plus_2",
     "ensembl_gene_id",
     "gene_biotype",
     "external_gene_name"),
  filter = "affy_hg_u133_plus_2",
  values = ids, uniqueRows=TRUE)

results=aggregate(external_gene_name ~., annotLookup, toString)
name = na.omit(results[match(ids, results$affy_hg_u133_plus_2),])

data("ALL")
dat <- exprs(ALL)
affyids = rownames(dat)
affyids


ovary = data.frame(name = ids, z_score= z_score,s=s, barcode=barcode)
ovary

affyid = getBM(c('affy_hg_u133_plus_2','external_gene_name'), filters = 'chromosome_name',values=list(1:23), mart=mart)
affyid
write.csv(affyid,"D:/Data Science/Data Science Core/Documents/Data Mining/Project/Ovary/affyid.csv", row.names = FALSE)



# 7. Export data to file.csv

write.csv(ovary,"D:/Data Science/Data Science Core/Documents/Data Mining/Project/Ovary/ovary.csv", row.names = FALSE)
write.csv(name,"D:/Data Science/Data Science Core/Documents/Data Mining/Project/Ovary/name.csv", row.names = FALSE)
gname <- read.csv(file="D:/Data Science/Data Science Core/Documents/Data Mining/Project/Gen name.csv", header=TRUE, sep=",")

