library(pheatmap)
library(PoiClaClu) # poissson Distance  
library(DESeq2)
library(ggplot2)
library(dplyr)
library(ggbeeswarm)
library(apeglm)
library(genefilter)
library(AnnotationDbi)
library(org.Hs.eg.db)

setwd("~/data_drive/Rnaseq_data")

res <- readRDS("EwS.rds")
counts <- assay(res)

meta <- as.data.frame (colData(res))


condition_df <- subset(meta,select = "condition")


# EDA analysis 

#Deseq 

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = condition_df,
                              design = ~condition)

dds <- DESeq(dds)
res <- results(dds,contrast = c("condition","shEF1","shCTR"))
res

summary(res)


# get the signnificant gene for down-regulation
resSig <- subset(res,padj <0.1)
head(resSig[ order(resSig$log2FoldChange), ])

# get the signicant gene for up -regulation 
head(resSig[ order(resSig$log2FoldChange, decreasing = TRUE), ])


# get the top gene  and plot using the plot counts 

topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, gene = topGene, intgroup=c("condition"))

# get the top gene and plot using the ggplot 
geneCounts <- plotCounts(dds, gene = topGene, intgroup = c("condition"),
                         returnData = TRUE)
ggplot(geneCounts, aes(x = condition, y = count, color = condition)) +
  scale_y_log10() +  geom_beeswarm(cex = 2)

ggplot(geneCounts, aes(x = condition, y = count, color = condition, group = condition)) +
  scale_y_log10() + geom_point(size = 3) + geom_line()


# MA plot by three techniques 
resultsNames(dds)
res <- lfcShrink(dds, coef="condition_shEF1_vs_shCTR", type="apeglm")
plotMA(res, ylim = c(-5, 5))

resnos <- lfcShrink(dds, coef="condition_shEF1_vs_shCTR")
plotMA(resnos, ylim = c(-5, 5))

plotMA(res, ylim = c(-5,5))
topGene <- rownames(res)[which.min(res$padj)]
with(res[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})

# histogram of p values
hist(res$pvalue[res$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")

#Gene Clustering
# calculating the variance stabilizing transfromation 
vsd  <- vst(dds,blind = FALSE)
head(assay(vsd),3)

topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 20)
mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("condition")])
rownames(anno) <- colnames(vsd)
colnames(anno)[1]<- "Condition"
pheatmap(mat, annotation_col = anno)


# Annotation using org.Hs.eg.db
columns(org.Hs.eg.db)
ens.str <- substr(rownames(res), 1, 15)
res$symbol <- mapIds(org.Hs.eg.db,
                     keys=ens.str,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
res$entrez <- mapIds(org.Hs.eg.db,
                     keys=ens.str,
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")


resOrdered <- res[order(res$pvalue),]
head(resOrdered)

results_ <- as.data.frame(resOrdered)

