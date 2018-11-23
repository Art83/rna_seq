
library(DESeq2)

# reading the data
counts <- read.delim("simple_counts.txt")
rownames(counts) <- counts$Geneid
counts <- counts[,-1]
colnames(counts) <- c("zero_1", "zero_2", "half_hour_1", "half_hour_2")

# dataframe with explicit names of groups
samples <- data.frame(timepoint = c(rep("zero_min",2), rep("half_hour",2)))

# Deseq
ds <- DESeqDataSetFromMatrix(countData=counts, colData=samples, design=~timepoint)
colnames(ds) <- colnames(counts)

ds <- DESeq(ds)
res <- results(ds, c("timepoint","zero_min", "half_hour"))

# getting rid of genes without p-values
sum( is.na(res$pvalue) )
res <- res[ ! is.na(res$pvalue), ]
sig <- res[ which(res$padj < 0.01), ]
sig <- sig[ order(sig$padj), ]
sig <- as.data.frame(sig)
# 6039 genes, 2654 pval< 0.01
upregulated <- sig[sig$log2FoldChange > 0,]
downregulated <- sig[sig$log2FoldChange < 0,]

#1284 downregulated genes vs 1370 upregulated genes


## Apply regularized-log transform to counts
rld <- rlogTransformation(ds)

pdf('Plot PCA')
## Principal component analysis
plotPCA(rld, intgroup="timepoint")
dev.off()

## Heatmap of sample distances
library("gplots")   # If this fails, run: install.packages("gplots")
library("RColorBrewer")
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix( sampleDists )
colours <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
heatmap.2(sampleDistMatrix,trace="none", col=colours)

## Heatmap of  50 most variable genes
library("genefilter")
topVarGenes <- head(order(rowVars(assay(rld)), decreasing=TRUE), 50)

pdf("heatmap")
heatmap.2(assay(rld)[rownames(counts), ], scale="row",
          trace="none", dendrogram="column", margins=c(10, 10),
          col=colorRampPalette(rev(brewer.pal(9, "RdBu")))(255))
dev.off()

pdf("volcano_plot")
#Adjusted P values (FDR Q values)
res <- as.data.frame(res)
with(res, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", cex=1.0, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~Q~value)))

with(subset(res, padj<0.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(padj), pch=20, col="red", cex=0.5))

#Add lines for absolute FC>2 and P-value cut-off at FDR Q<0.01
abline(v=0, col="black", lty=3, lwd=1.0)
abline(v=-2, col="black", lty=4, lwd=2.0)
abline(v=2, col="black", lty=4, lwd=2.0)
abline(h=-log10(max(res$pvalue[res$padj<0.01], na.rm=TRUE)), col="black", lty=4, lwd=2.0)
dev.off()

write.table(rownames(sig[1:50,]), "50_most_expressed.txt", quote = F, row.names = F, col.names = F)
write.table(rownames(upregulated[1:50,]), "50_most_expressed_upregulated.txt", quote = F, row.names = F, col.names = F)
write.table(rownames(downregulated[1:50,]), "50_most_expressed_downregulated.txt", quote = F, row.names = F, col.names = F)
