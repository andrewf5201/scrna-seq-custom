library('properties') 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager") 

BiocManager::install("DESeq2")

library("DESeq2")
library('ggplot2') 
library('geneplotter')
 
library('vsn')

#get input values
args <- commandArgs(trailingOnly = TRUE)
algorithm <-args[1]
config_file <- args[2]

algorithm = "default"
config_file = "/Users/lwang1/pipeline/config/config_bc.properties"

#load configuration file
config <-read.properties (config_file)
genomeAnnotation= config$genome_type
pipeline.dir = config$pipeline_dir
experiment = config$experiment
gene.length.file = config$gene_length_file
count.matrix.file = config$matrix_file

#construct input and output directories
analysis.dir = paste(pipeline.dir,"analysis", experiment, sep="/")
input.folder = paste(analysis.dir, "data", sep="/")
output.folder = paste(analysis.dir, "results", "Deseq",sep="/")

 

# read input count files
setwd(input.folder)
#filelist = list.files(pattern="*counts.txt")
#print(filelist)

#assuming tab separated values with a header
#datalist = lapply(filelist, function(x)read.table( x, header=TRUE, row.names=1, fill = FALSE))
#cts = do.call("cbind", datalist)
#remove last 5 rows which contains summary info
#cts <- cts[1:(nrow(cts)-5),]

# load count matrix file
cts <- read.table(count.matrix.file, sep="\t", header=TRUE)
print(cts)

#only keep 2000 genes based on variance
library(countSubsetNorm)
cts=var_filter_counts(cts, n_keep=2000, genes="rows",log2_transform=FALSE)
print("--Genes trimmed by Variance")

df <- read.table(metadata.file, header = TRUE, sep = "\t")
#get all factors defined in metadata file (from the 6th column)
fctrs=colnames(df)[6:length(colnames(df))]
factorlists=c()
library(sqldf)

for (fctr in fctrs){
  print(fctr)
  column_list=c()
  for (i in (1:ncol(cts))){
    id=colnames(cts)[i]
    sourcedf <- fn$sqldf("select $fctr from df where ID='$id' ")
    sourcename <- as.character((sourcedf[1,1]))
    column_name=sourcename
    column_list=c(column_list,column_name)
  }
  factorlists[[length(factorlists)+1]]=column_list
}

#colData=cbind(gsub("\\_.*","",colnames(cts)), conditionList )
colData=as.data.frame(factorlists)
colnames(colData)=fctrs
rownames(colData)=colnames(cts)
print("--colData created")

designvar=paste(fctrs, collapse =" + ")
designvar=as.formula(paste('~',designvar,sep=""))

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = colData,
                              design= designvar)

#dds <- DESeqDataSetFromMatrix(countData = d,
#                              colData = ann,
#                            design= ~cell_type1)
print("--DESeqDataSet created")
keep <- rowSums(counts(dds)) >= 25
dds <- dds[keep,]

vsd <- vst(dds, blind=FALSE)
print("--vst completed")
head(assay(vsd), 3)


sdPlot = paste(output.folder, '/meanSd.png', sep="")
png(sdPlot)
meanSdPlot(assay(vsd))
dev.off()

if (algorithm == "default") {
  print("--default algorithm")
  dds <- estimateSizeFactors(dds)
  head(sizeFactors(dds))
} else if (algorithm == "tmm") {
  print("--TMM")
  library(edgeR)
  par(mfrow=c(1,1))
  dgList <- DGEList(counts=cts, genes=rownames(cts))
  dgList <- calcNormFactors(dgList, method="TMM")
  plotMDS.default(dgList)
  #sizeFactors(dds)=dgList$samples$lib.size * dgList$samples$norm.factors
  sizeFactors(dds)=dgList$samples$norm.factors
} else if (algorithm == "fpkm") {
  print("--FPKM")
  #gene_len <-read.delim(paste(input.folder,"hg38_gene_length.txt",sep=""), header=F, sep="\t")
  gene_len <- read.delim(gene.length.file, header=F, sep="\t")
  colnames(gene_len) <- c("GeneName","Len")
  mcols(dds)$basepairs <- gene_len$Len
  fpkm(dds, gene.lengths)
}
print("--Normalization complete")

normPlot1 = paste(output.folder, '/normalization.png', sep="")
png(normPlot1)
par(mfrow=c(1,1))
plot( log2( 1 + counts(dds, normalized=TRUE)[ , 1:2] ),col=rgb(0,0,0,.2), pch=16, cex=0.3 )
dev.off()

normPlot2 = paste(output.folder, '/normalization_raw.png', sep="")
png(normPlot2)
plot( log2( 1 + counts(dds)[ , 1:2] ), col=rgb(0,0,0,.2), pch=16, cex=0.3 )
dev.off()

normPlot3 = paste(output.folder, '/normalization_vsd.png', sep="")
png(normPlot3)
plot( assay(vsd)[ , 1:2], col=rgb(0,0,0,.2), pch=16, cex=0.3 )
dev.off()
print("-- Finish geneating normalization plot")


notAllZero <- (rowSums(counts(dds))>0)
# Non-Zero
random=base::sample(x=ncol(counts(dds)),size=25,replace=FALSE)

meanSd1 = paste(output.folder, '/meanSd_norm.png', sep="")
png(meanSd1)
meanSdPlot(log2(counts(dds,normalized=TRUE)[notAllZero,] + 1))
dev.off()

meanSd2 = paste(output.folder, '/meanSd_raw.png', sep="")
png(meanSd2)
meanSdPlot(log2(counts(dds)[notAllZero,] + 1))
dev.off()

meanSd3 = paste(output.folder, '/meanSd_vsd.png', sep="")
png(meanSd3)
meanSdPlot(assay(vsd[notAllZero,]))
dev.off()

boxplot = paste(output.folder, '/boxplot_raw.png', sep="")
png(boxplot)
boxplot(log2(counts(dds)[notAllZero,random] + 1))
dev.off()

boxplot2 = paste(output.folder, '/boxplot_norm.png', sep="")
png(boxplot2)
boxplot(log2(counts(dds,normalized=TRUE)[notAllZero,random] + 1))
dev.off()

print("-- Finish geneating meanSd plot")

library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=FALSE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,fctrs])
rownames(df)<-colnames(assay(vsd)[select,])
heatmapp = paste(output.folder, '/pheatmap.png', sep="")
png(heatmapp)
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
dev.off()
print("--pheatmap created")


sampleDists <- dist(t(assay(vsd)))
plot(hclust(sampleDists,method="ward.D"))


library(RColorBrewer)

sampleDistMatrix <- as.matrix( sampleDists )
#rowname=paste("vsd$",fctrs,sep="")
#rowlist=c()
#for (i in rowname){
#  rowlist=c(rowlist,parse(text=i))
#}
#rowlist=eval(rowlist)
#rownames(sampleDistMatrix) <- paste(rowlist, sep="-" )
rownames(sampleDistMatrix) <- paste(vsd$molecular_subtype, sep="-" )
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
library("rafalib")
cols=palette(brewer.pal(8, "Dark2"))[as.fumeric(as.character(colData[specifiedFactor]))]
hc <- hclust(sampleDists, method="ward.D")
library(gplots)
heatmap = paste(output.folder, '/heatmap.png', sep="")
png(heatmap)
heatmap.2( sampleDistMatrix, Rowv=as.dendrogram(hc),
           symm=TRUE, trace="none", col=colors,
           ColSideColors = cols,
           margins=c(1,5), labCol=FALSE )
dev.off()
library(dendextend)

dend <- as.dendrogram(hc)
# Like:
# dend <- USArrests[1:5,] %>% dist %>% hclust %>% as.dendrogram

dend1 <- color_branches(dend, k = 4)
colors_to_use <- as.numeric(colData[specifiedFactor])
colors_to_use <- colors_to_use[order.dendrogram(dend)]
labels_colors(dend1) <- colors_to_use


plot(dend1, main = "Colored branches")

# install.packages("dendextend")
library(dendextend)

dend <- as.dendrogram(hc)
# Like:
# dend <- USArrests[1:5,] %>% dist %>% hclust %>% as.dendrogram

dend1 <- color_branches(dend, k = 4)
colors_to_use <- as.numeric(colData[specifiedFactor])
colors_to_use <- colors_to_use[order.dendrogram(dend)]
labels_colors(dend1) <- colors_to_use


plot(dend1, main = "Colored branches")




library(psych)
sampleSpearman=cor(assay(vsd),method="spearman")
sampleSpearmanDist=cor2dist(sampleSpearman)

rownames(sampleSpearmanDist) <- paste(vsd$molecular_subtype, sep="-" )
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
plot(hclust(as.dist(sampleSpearmanDist), method="ward.D2"))
cols=palette(brewer.pal(8, "Dark2"))[as.fumeric(as.character(colData[specifiedFactor]))]
hc <- hclust(as.dist(sampleSpearmanDist), method="ward.D")
library(gplots)
heatmap = paste(output.folder, '/heatmapSpearman.png', sep="")
png(heatmap)
heatmap.2(sampleSpearmanDist, Rowv=as.dendrogram(hc),
          symm=TRUE, trace="none", col=colors,
          ColSideColors = cols,
          margins=c(1,5), labCol=FALSE )
dev.off()
library(dendextend)



dend <- as.dendrogram(hc)
# Like:
# dend <- USArrests[1:5,] %>% dist %>% hclust %>% as.dendrogram

dend1 <- color_branches(dend, k = 4)
colors_to_use <- as.numeric(colData[specifiedFactor])
colors_to_use <- colors_to_use[order.dendrogram(dend)]
labels_colors(dend1) <- colors_to_use


plot(dend1, main = "Colored branches")


samplePearson=cor(assay(vsd),method="pearson")
samplePearsonDist=cor2dist(samplePearson)

rownames(samplePearsonDist) <- paste(vsd$molecular_subtype, sep="-" )
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
plot(hclust(as.dist(samplePearsonDist)))
cols=palette(brewer.pal(8, "Dark2"))[as.fumeric(as.character(colData[specifiedFactor]))]
hc <- hclust(as.dist(samplePearsonDist), method="ward.D")
plot(hc)
clusters=cutree(hc,k=4)

library(gplots)
heatmap = paste(output.folder, '/heatmap.png', sep="")
png(heatmap)
heatmap.2(samplePearsonDist, Rowv=as.dendrogram(hc),
          symm=TRUE, trace="none", col=colors,
          ColSideColors = cols,
          margins=c(1,5), labCol=FALSE )
dev.off()

print("--Heatmaps created")

PCA = paste(output.folder, '/PCA.png', sep="")
png(PCA, width=800, height=800)
par(mfrow=c(1,1))
pca <-plotPCA(vsd, intgroup=fctrs, returnData=TRUE)
percentVar <- round(100 * attr(pca, "percentVar"))
if(length(fctrs)==2){
  ggplot(pca, aes(PC1, PC2, color=fctrs[1], shape=fctrs[2])) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    coord_fixed(ratio=3 )+
    theme(legend.text=element_text(size=12))
}else {
  ggplot(pca, aes(PC1, PC2, color=molecular_subtype)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    coord_fixed( ratio=3 )+
    theme(legend.text=element_text(size=10))

}
dev.off()


plotPCA(vsd, intgroup=fctrs, ntop=100)

print("--PCA Completed")

dds <- DESeq(dds, test="LRT", reduced=~1, sfType="poscounts", useT=TRUE, minmu=1e-6, minReplicatesForReplace=Inf)
res <- results(dds)

print("DESeq completed")

#library("IHW")
#resIHW <- results(dds, filterFun=ihw)
#summary(resIHW)

resultsNames(dds)

resLFC <- lfcShrink(dds, coef=resultsNames(dds)[2], type="apeglm")
resLFC

summary(res)

resOrdered <- res[order(res$pvalue),]
head(resOrdered)
write.table(resOrdered,"resOrdered.txt",row.names = TRUE, col.names = TRUE)

random=base::sample(x=ncol(assays(dds)[["cooks"]]),size=10,replace=FALSE)
boxplot(log10(assays(dds)[["cooks"]])[,random], range=0, las=2)



#library("Rfast")
#cvmatrix=rowcvs(counts(dds))
#cvgenematrix=cbind(rownames(counts(dds)),cvmatrix)
#cvOrdered=cvgenematrix[order(cvgenematrix[,2],decreasing = TRUE)]

sampleDists <- dist(t(assay(vsd[rownames(resOrdered)[1:100],])))
plot(hclust(sampleDists, method="ward.D"))
sampleDistMatrix <- as.matrix( sampleDists )
#rowname=paste("vsd$",fctrs,sep="")
#rowlist=c()
#for (i in rowname){
#  rowlist=c(rowlist,parse(text=i))
#}
#rowlist=eval(rowlist)
#rownames(sampleDistMatrix) <- paste(rowlist, sep="-" )
rownames(sampleDistMatrix) <- paste(vsd$molecular_subtype, sep="-" )
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
library("rafalib")
cols=palette(brewer.pal(8, "Dark2"))[as.fumeric(as.character(colData[specifiedFactor]))]
hc <- hclust(sampleDists, method="ward.D")
library(gplots)
heatmap = paste(output.folder, '/heatmap.png', sep="")
png(heatmap)
heatmap.2( sampleDistMatrix, Rowv=as.dendrogram(hc),
           symm=TRUE, trace="none", col=colors,
           ColSideColors = cols,
           margins=c(1,5), labCol=FALSE )
dev.off()
library(dendextend)

dend <- as.dendrogram(hc)
# Like:
# dend <- USArrests[1:5,] %>% dist %>% hclust %>% as.dendrogram

dend1 <- color_branches(dend, k = 4)
colors_to_use <- as.numeric(colData[specifiedFactor])
colors_to_use <- colors_to_use[order.dendrogram(dend)]
labels_colors(dend1) <- colors_to_use

library(psych)
sampleSpearman=cor(assay(vsd)[rownames(resOrdered)[1:1000],],method="spearman")
sampleSpearmanDist=cor2dist(sampleSpearman)

rownames(sampleSpearmanDist) <- paste(vsd$molecular_subtype, sep="-" )
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
plot(hclust(as.dist(sampleSpearmanDist),method="ward.D"))
cols=palette(brewer.pal(8, "Dark2"))[as.fumeric(as.character(colData[specifiedFactor]))]
hc <- hclust(as.dist(sampleSpearmanDist), method="ward.D")
library(gplots)
heatmap = paste(output.folder, '/heatmap.png', sep="")
png(heatmap)
heatmap.2(sampleSpearmanDist, Rowv=as.dendrogram(hc),
          symm=TRUE, trace="none", col=colors,
          ColSideColors = cols,
          margins=c(1,5), labCol=FALSE )

library(dendextend)



dend <- as.dendrogram(hc)
# Like:
# dend <- USArrests[1:5,] %>% dist %>% hclust %>% as.dendrogram

dend1 <- color_branches(dend, k = 4)
colors_to_use <- as.numeric(colData[specifiedFactor])
colors_to_use <- colors_to_use[order.dendrogram(dend)]
labels_colors(dend1) <- colors_to_use


plot(dend1, main = "Colored branches")


samplePearson=cor(assay(vsd)[rownames(resOrdered)[1:100],],method="pearson")
samplePearsonDist=cor2dist(samplePearson)

rownames(samplePearsonDist) <- paste(vsd$molecular_subtype, sep="-" )
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

cols=palette(brewer.pal(8, "Dark2"))[as.fumeric(as.character(colData[,1]))]
hc <- hclust(as.dist(samplePearsonDist), method="ward.D")
plot(hc)
clusters=cutree(hc,k=4)

library(gplots)
heatmap = paste(output.folder, '/heatmap.png', sep="")
png(heatmap)
heatmap.2(samplePearsonDist, Rowv=as.dendrogram(hc),
          symm=TRUE, trace="none", col=colors,
          ColSideColors = cols,
          margins=c(1,5), labCol=FALSE )
dev.off()
# install.packages("dendextend")
library(dendextend)

dend <- as.dendrogram(hc)
# Like:
# dend <- USArrests[1:5,] %>% dist %>% hclust %>% as.dendrogram

dend1 <- color_branches(dend, k = 4)
colors_to_use <- as.numeric(colData[specifiedFactor])
colors_to_use <- colors_to_use[order.dendrogram(dend)]
labels_colors(dend1) <- colors_to_use


plot(dend1, main = "Colored branches")






print("--Heatmaps created")

PCAres = paste(output.folder, '/PCAres.png', sep="")
png(PCAres, width=800, height=800)
par(mfrow=c(1,1))
pca <-plotPCA(vsd[rownames(resOrdered)[1:100]], intgroup=fctrs, returnData=TRUE)
percentVar <- round(100 * attr(pca, "percentVar"))
if(length(fctrs)==2){
  ggplot(pca, aes(PC1, PC2, color=fctrs[1], shape=fctrs[2])) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    coord_fixed(ratio=3 )+
    theme(legend.text=element_text(size=12))
}else {
  ggplot(pca, aes(PC1, PC2, color=molecular_subtype)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    coord_fixed( ratio=3 )+
    theme(legend.text=element_text(size=10))

}
dev.off()


plotDispEsts(dds)

plotMA(resLFC, ylim=c(-10,20))

topGene <- rownames(res)[which.min(res$padj)]
par( mfrow = c( 1, 1 ) )
plotCounts(dds, gene=topGene, intgroup=fctrs)
DESeq2::plotMA(res, ylim=c(-10,20))
with(res[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})
plotDispEsts(dds)
hist(res$pvalue, breaks=20, col="grey50", border="white")
hist(res$pvalue[res$baseMean > 1], breaks=20, col="grey50", border="white")

library("genefilter")
topVarGenes <- head(order(-rowVars(assay(vsd))),35)

colors <- colorRampPalette( rev(brewer.pal(9, "PuOr")) )(255)
sidecols <- c("grey","dodgerblue")[ vsd$molecular_subtype ]
mat <- assay(vsd)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
heatmapTVG = paste(output.folder, '/heatmapTVG.png', sep="")
png(heatmapTVG)
heatmap.2(mat, trace="none", col=colors, ColSideColors=sidecols,
          labRow=FALSE, mar=c(10,2), scale="row")
dev.off()

library(ggplot2)
library(Rtsne)
tsne=Rtsne(t(assay(vsd)), perplexity=5)

tsneplot = paste(output.folder, '/tsne5.png', sep="")
png(tsneplot)
tsne_plot <- data.frame(x = tsne$Y[,1], y = tsne$Y[,2], col = vsd$molecular_subtype)
ggplot(tsne_plot) + geom_point(aes(x=x, y=y, color=col))+ggtitle("Perplexity = 5")
dev.off()

tsneRes=Rtsne(t(assay(vsd[rownames(resOrdered)[1:100]])), perplexity=5)

tsneResPlot = paste(output.folder, '/tsne5resOrdered.png', sep="")
png(tsneResPlot)
tsne_plot <- data.frame(x = tsneRes$Y[,1], y = tsneRes$Y[,2], col = vsd$molecular_subtype)
ggplot(tsne_plot) + geom_point(aes(x=x, y=y, color=col))+ggtitle("Perplexity = 5")
dev.off()


library("AnnotationDbi")
library("org.Hs.eg.db")
library("org.Mm.eg.db")
convertIDs <- function( ids, from, to, db, ifMultiple=c("putNA", "useFirst")) {
  stopifnot( inherits( db, "AnnotationDb" ) )
  ifMultiple <- match.arg( ifMultiple )
  suppressWarnings( selRes <- AnnotationDbi::select(
    db, keys=ids, keytype=from, columns=c(from,to) ) )
  if ( ifMultiple == "putNA" ) {
    duplicatedIds <- selRes[ duplicated( selRes[,1] ), 1 ]
    selRes <- selRes[ ! selRes[,1] %in% duplicatedIds, ]
  }
  return( selRes[ match( ids, selRes[,1] ), 2 ] )
}
tmp=gsub("\\..*","",row.names(res))
if ( genome=="human" ) {
  genomeAnnotation= org.Hs.eg.db
}else if ( genome=="mouse" ){
  genomeAnnotation= org.Mm.eg.db
}

res$symbol <- mapIds(genomeAnnotation,
                     keys=tmp,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
print("--Gene Symbols added")

res$entrez <- mapIds(genomeAnnotation,
                     keys=tmp,
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")
print("--ENTREZID ID's added")

resOrdered <- res[order(res$pvalue),]
head(resOrdered)
resOrdered=na.omit(resOrdered)
write.table(resOrdered,"resOrdered.txt",row.names = TRUE, col.names = TRUE)
na
