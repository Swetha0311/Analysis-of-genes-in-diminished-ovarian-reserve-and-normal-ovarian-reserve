library (DESeq2)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)

setwd("/Users/swethayadavalli/Desktop/wrk")
count_new_data=read.table("dor_nor.csv", sep=",", header=T)

class(count_new_data)
head(count_new_data)

nrow (count_new_data)
sum(duplicated(count_new_data))

print("Row and col positions of NA values")
which(is.na(count_new_data))

genenames <- count_new_data$gene_name
head(genenames)

count_new_data <- count_new_data [,2:7]
rownames(count_new_data) <- genenames

count_new_data <- as.matrix(count_new_data)
class(count_new_data)

coldata <- data.frame("condition"=as.factor (c(rep ("Diminished",3),rep ("Normal",3))), row.names=colnames (count_new_data))


count_new_data<-count_new_data [rowSums (count_new_data) >= 5, ]
head(count_new_data)

dds <- DESeqDataSetFromMatrix(countData = count_new_data, colData = coldata, design= ~ condition)
dds

dds <- DESeq(dds)

normalized_counts <- counts (dds, normalized = TRUE)
head (normalized_counts)
write.csv (normalized_counts, 'normalized_counts.csv')

res <- results(dds, contrast = c('condition', 'Diminished', 'Normal'), alpha = 0.05)
res

summary(res)


resordered <- res [order(res$padj),]
resordered

plotMA (res,cex =0.7, ylim=c(-20,20))
abline(h=c(-1,1), col="red", lwd=3)

resultsNames (dds)
resLFC <- lfcShrink(dds, coef="condition_Normal_vs_Diminished", type = "apeglm")
plotMA (resLFC,cex =0.7, ylim=c(-10,10))
abline (h=c(-1,1), col="red", lwd=3)


plotDispEsts(dds, main = "Dispersion Plot")

rld <- rlogTransformation (dds, blind = FALSE)
head (assay(rld))
hist (assay(rld))


PCAA <- plotPCA(rld, intgroup='condition')
PCAA + geom_text (aes (label = name),size=2.5)+ggtitle('PCA Plot')


EnhancedVolcano (res,
  lab=rownames (res),
  x ='log2FoldChange',
  y = 'pvalue')
EnhancedVolcano (res,
  lab = rownames(res),
  x = 'log2FoldChange',
  y = 'pvalue',
  title = 'Normal vs. Diminished',
  pCutoff = 10e-4,
  FCcutoff = 0.5,
  pointSize = 3.0,
  labSize = 6.0)


sampleDists <- dist(t(assay(rld)))
library(RColorBrewer)
sampleDistMatrix <- as.matrix (sampleDists)
colnames (sampleDistMatrix)

colors <- colorRampPalette(rev(brewer.pal (9, "Blues")) )(255)
pheatmap (sampleDistMatrix,
          clustering_distance_rows=sampleDists,
          clustering_distance_cols=sampleDists,
          col=colors)

top_genes <- res [order (res$padj),] [1:20,]
class (top_genes)
top_genes
top_genes <- row.names(top_genes)
top_genes
pheatmap (assay(rld) [top_genes,],cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=TRUE, fontsize_row =8,
          annotation_col = coldata)
write.csv(top_genes, file = "top_genes.csv")