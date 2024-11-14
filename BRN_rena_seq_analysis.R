library(pheatmap)
library(PoiClaClu) # poissson Distance  
library(DESeq2)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(ggbeeswarm)
library(apeglm)
library(genefilter)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(EnhancedVolcano)
library(msigdbr)
library(clusterProfiler)
library(tibble)

setwd("~/data_drive/Rnaseq_data")

res <- readRDS("EwS.rds")
counts <- assay(res)

meta <- as.data.frame (colData(res))


condition_df <- subset(meta,select = "condition")



#Deseq 

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = condition_df,
                              design = ~condition)

nrow(dds)


vsd  <- vst(dds,blind = FALSE)
head(assay(vsd),3)

# PCA Plot to view the sample variance within the dataset 
# PCA plot using ggplot 

pcaData <- plotPCA(vsd, intgroup = c( "condition"), returnData = TRUE)
pcaData
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = condition, shape = condition)) +
  geom_point(size =2) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data")+
  theme(plot.title = element_text(hjust = 0.5))
  
dds <- DESeq(dds)
res <- results(dds,contrast = c("condition","shEF1","shCTR"))
res

summary(res)



# Mean count and log fold change 
resultsNames(dds)
res <- lfcShrink(dds, coef="condition_shEF1_vs_shCTR", type="apeglm")
plotMA(res, ylim = c(-5, 5))


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
res



EnhancedVolcano(res,
                lab = res$symbol,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'shCTR versus shEF1')

# getting the significiant genes 
resSig <- res[which (res$padj < 0.01 & abs(res$log2FoldChange) >= 1 & res$baseMean >= 20),]
resSig$Ensemble <- substr(rownames(resSig), 1, 15)
rld <-  rlog(dds)
mat <- assay(rld)



# get the DE genes with Downregulation 
orderSig <- resSig[order(resSig$padj),]
id <- rownames(orderSig)
id1 <- orderSig$Ensemble
id2 <- orderSig$symbol
DE <- mat[id,]
rownames(DE) <- id1
rownames(DE)<- id2
top10DE <- head(DE,n=10)
pheatmap(top10DE,scale = "row",clustering_distance_rows = "correlation",annotation_col = condition_df, main = "Top 10 Differently Expresssed Genes")


# Get top 10 downregulated genes (by most significant p-value)
downregulated <- orderSig[orderSig$log2FoldChange < 0, ]
down10DE <- head(downregulated, n=10)
down10DE_mat <- mat[rownames(down10DE), ]
rownames(down10DE_mat) <- down10DE$symbol

# Plot heatmap for downregulated genes
pheatmap(down10DE_mat, scale="row", clustering_distance_rows="correlation", 
         annotation_col=condition_df, main="Top 10 Downregulated Genes")



# --------------------------------------------- Enrichemnt analysis  -------------------------------------------
genes_set <- msigdbr(species = "Homo sapiens",category = "C5")
genes_set <- genes_set %>%
  dplyr::select(gs_name,gene_symbol)


overexpressed_genes <- as.data.frame(res) %>%
  dplyr::filter(padj < 0.1 & log2FoldChange > 2) %>%
  pull(symbol)


over_expressed <- enricher(gene=overexpressed_genes,
                           TERM2GENE = genes_set)
a <- as.data.frame(over_expressed)
dotplot(over_expressed)
barplot(over_expressed)


downregulated_genes <- as.data.frame(res) %>%
  dplyr::filter(padj < 0.1 & log2FoldChange < 0) %>%
  pull(symbol)

downregulated_expressed <- enricher(gene=downregulated_genes,
                           TERM2GENE = genes_set)
b <- as.data.frame(downregulated_expressed)
dotplot(downregulated_expressed)
barplot(downregulated_expressed)





#Enrichemnt analysis 
resultdf <- as.data.frame(res) %>%
  arrange(padj) %>%  # Sorts the data by the padj column in ascending order
  mutate(gse_metric = -log10(padj) * sign(log2FoldChange))  # Calculates gse_metric

# dealing with the infinity values 
resultdf <- resultdf  %>%
  mutate(padj=case_when(padj==0 ~ .Machine$double.xmin,
                        TRUE ~ padj)) %>%
  mutate(gse_metric= -log10(padj) * sign(log2FoldChange))



#GSea --> histogram
hist(resultdf$gse_metric,breaks = 100)

ranks <- resultdf %>%
  filter(!is.na(gse_metric))%>%
  filter(!is.na(symbol))%>%
  arrange(desc(gse_metric)) %>%
  dplyr::select(symbol,gse_metric) %>%
  distinct(symbol, .keep_all = TRUE) %>%
  deframe()

gseares <- GSEA(geneList = ranks,TERM2GENE = genes_set)

GSEA_df <- as.data.frame(gseares)

dotplot(gseares)


# get the top over expressed pathways 
top_pathways <- GSEA_df %>%
  top_n(n=10,wt = NES) %>%
  pull(ID)

top__pathways_plots <- lapply(top_pathways, function(pathways){
  gseaplot(gseares,geneSetID = pathways,title = pathways)
})

print(top__pathways_plots)

# get the bottom expressed pathways 

bottom_pathways <- GSEA_df %>%
  top_n(n=10,wt = -NES) %>%
  pull(ID)

bottom_pathways_plot <- lapply(bottom_pathways, function(pathways){
  gseaplot(gseares,geneSetID = pathways,title = pathways)
})

print(bottom_pathways_plot)


