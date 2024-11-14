library(pheatmap)
library(PoiClaClu) # poissson Distance  
library(DESeq2)
library(ggplot2)
library(dplyr)

res <- readRDS("EwS.rds")
counts <- assay(res)
 
meta <- as.data.frame (colData(res))


condition_df <- subset(meta,select = "condition")


# EDA analysis 

#Deseq 

dds <- DESeqDataSetFromMatrix(countData = counts,
                             colData = condition_df,
                             design = ~condition)

nrow(dds)

# filtering the dataset 
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
nrow(dds)


# calculating the variance stabilizing transfromation 
vsd  <- vst(dds,blind = FALSE)
head(assay(vsd),3)

# calculating the rlog
rld <- rlog(dds,blind=FALSE)
head(assay(rld),3)

dds <- estimateSizeFactors(dds)
df <- bind_rows(
  as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))

colnames(df)[1:2] <- c("x", "y")  

lvls <- c("log2(x + 1)", "vst", "rlog")
df$transformation <- factor(df$transformation, levels=lvls)

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)



sampleDists <- dist(t(assay(vsd)))
sampleDists

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( dds$condition)
colnames(sampleDistMatrix) <- NULL
colors <- colorRamp(c("red", "blue"))
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         main="Clustering  Samples Based VST ")



# Posssion Distance based sample clustering 
poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix( poisd$dd)
rownames(samplePoisDistMatrix) <- paste( dds$condition, sep = " - " )
colors <- colorRamp(c("red", "blue"))
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         main="Clustering  Samples Based on Possion Distance")



# PCA plot using deseq
DESeq2::plotPCA(vsd,intgroup = c("condition"))
#PCA plot using ggplot 
pcaData <- plotPCA(vsd, intgroup = c( "condition"), returnData = TRUE)
pcaData
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(x = PC1, y = PC2, color = condition, shape = condition)) +
  geom_point(size =2) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data")


#MDS plot 

mds <- as.data.frame(colData(vsd))  %>%
  cbind(cmdscale(sampleDistMatrix))
ggplot(mds, aes(x = `1`, y = `2`, color = condition, shape = condition)) +
  geom_point(size = 2) + coord_fixed() + ggtitle("MDS with VST data")

# 
mdsPois <- as.data.frame(colData(dds)) %>%
  cbind(cmdscale(samplePoisDistMatrix))
ggplot(mdsPois, aes(x = `1`, y = `2`, color = condition, shape = condition)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS with PoissonDistances")

