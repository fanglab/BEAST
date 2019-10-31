#!/usr/bin/env Rscript

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")

library("DESeq2")
library("RColorBrewer")
library("gplots")
library("ggplot2")
library("rgl")
library("pheatmap")

args <- commandArgs()

#Run DESeq2
countData <- as.matrix(read.csv(args[6],header=TRUE,sep="\t",row.names="gene_id"))
#countData <- as.matrix(read.csv("count_file",header=TRUE,sep="\t",row.names="gene_id"))
countData <- countData[ rowSums(countData>1)>= 6,] #Select genes with at least one count in all samples
colData <- read.csv(args[7], row.names="file")
#colData <- read.csv("colData.csv", row.names="file")
condition <- factor(c(rep("untreated",3),rep("treated",3)), levels=c("untreated","treated"))
colData$condition <- relevel(colData$condition, "untreated") #Here we are telling that we want the control sample to be the 'untreated' one.

dds <- DESeqDataSetFromMatrix(countData = countData,colData = colData,design = ~condition)
dds <- dds[ rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds)
res <- results(dds)
reshrunk <- lfcShrink(dds, contrast=c("condition","treated","untreated"), res=res)

comparison <- paste( levels(dds$condition)[2] , 'VS' , levels(dds$condition)[1] , sep='_')
normalized.data <- counts(dds,normalized=TRUE)  # divided by sizeFactors(dds)
log2fc <- res$log2FoldChange #log2 fold change (MLE)
lfcSE	<- res$lfcSE #standard error
pvalue <- res$pvalue #Wald test p-value
fdr <- res$padj  #False discovery rate

outInfo <- data.frame( countData, normalized.data , log2fc , lfcSE , pvalue , fdr )
colnames(outInfo) = c( paste( 'original',condition,sep='.' ) , paste( 'normalized',condition,sep='.') , paste(comparison,'log2fc',sep='_'), 'lfcSE' , 'pvalue' , 'fdr' ) 
cat( sum(fdr<0.01,na.rm=T) , 'genes with fdr<0.01' )
write.csv(outInfo, file = "outInfo.txt")

plotMA(res, ylim=c(-2,2))
dev.copy(pdf,'plotMA.pdf')
dev.off()
plotMA(reshrunk, ylim=c(-2,2))
dev.copy(pdf,'plotMAShrunk.pdf')
dev.off()
write.csv(res, file = "MyData.csv")

#Make dendrogram to look at sample clustering and 2D PCA plot to check for potential outliers
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
rld<-rlog(dds, blind=F) #To ensure we have a roughly equal contribution from all genes, we use rlog-transformed data
distsRL <- dist(t(assay(rld))) #Euclidean distances between samples. We need to transpose the matrix of values using t, because the dist function expects the different samples to be rows of its argument, and different dimensions (here, genes) to be columns.
mat <- as.matrix(distsRL)
#rownames(mat) <- paste(condition,countData , sep=' : '))
hc <- hclust(distsRL)
heatmap.2(mat, Rowv=as.dendrogram(hc),
          symm=TRUE, trace='none',
          col = rev(hmcol), margin=c(13, 13))
dev.copy(pdf,'dendrogram.pdf')
dev.off()

print(plotPCA(rld, intgroup=c("condition")))
dev.copy(pdf,'PCA.pdf')
dev.off()

