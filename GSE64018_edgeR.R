# -------------------------- RNA-Seq of ASD vs control -------------------------
library(edgeR)
library(org.Hs.eg.db)
# --------------------------------- GSE64018 -----------------------------------
## 1. Reading in the data
rawdata_ASD <- read.delim("GSE64018_countlevel_12asd_12ctl.txt")
head(rawdata_ASD)
group <- factor(c(rep("ASD",12),rep("control",12)))

y <- DGEList(counts=rawdata_ASD, genes=row.names(rawdata_ASD), group = group)

## 2. Annotation

### Retain only those transcripts with IDs in the current NCBI annotation, 
### which is provided by the org.HS.eg.db package:
idfound <- y$genes$genes %in% mappedRkeys(org.Hs.egENSEMBL)
y <- y[idfound,]
dim(y)

### We add Entrez Gene IDs to the annotation:
egENSEMBL <- toTable(org.Hs.egENSEMBL)
head(egENSEMBL)
m <- match(y$genes$genes, egENSEMBL$ensembl_id)
y$genes$EntrezGene <- egENSEMBL$gene_id[m]

### Now use the Entrez Gene IDs to update the gene symbols:
egSYMBOL <- toTable(org.Hs.egSYMBOL)
head(egSYMBOL)
m <- match(y$genes$EntrezGene, egSYMBOL$gene_id)
y$genes$Symbol <- egSYMBOL$symbol[m]
head(y$genes)

## 3. Filtering and normalization

### Different RefSeq transcripts for the same gene symbol count predominantly the
### same reads. So we keep one transcript for each gene symbol. We choose the 
### transcript with highest overall count:

o <- order(rowSums(y$counts), decreasing=TRUE)
y <- y[o,]
d <- duplicated(y$genes$Symbol)
y <- y[!d,]
nrow(y)

### A gene is only retained if it is expressed at a minimum level:
keep <- filterByExpr(y)
summary(keep)
y <- y[keep, , keep.lib.sizes=FALSE]

### Apply TMM normalization
y <- calcNormFactors(y)
y$samples

## Data exploration
### Generate MDS plot
color_easy = c("red","blue")[group] # assign color to factors (ASD vs control)
plotMDS(y, col=color_easy, cex = 0.5)

design <- model.matrix(~group)
y <- estimateDisp(y,design)

## Perform quasi-likelihood F-tests:

fit <- glmQLFit(y,design)
qlf <- glmQLFTest(fit,coef=2)
save(qlf, file = 'GSE64018_ASD.qlf.RData')

summary(decideTests(qlf, adjust.method = 'none', p.value = 0.01))
ASD_deGenes <- topTags(qlf, n=1218+364)
ASD_deGenes <- ASD_deGenes$table
save(ASD_deGenes, file='GSE64018_ASD_deGenes.RData') 

ASDgenes <- topTags(qlf, n=19000)
ASDgenes <- ASDgenes$table
save(ASDgenes, file='GSE64018_ASD_allGenes.RData') 