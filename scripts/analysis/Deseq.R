install.packages("remotes")
remotes::install_github("mjdufort/countSubsetNorm")
library(countSubsetNorm)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager") 

install.packages("rafalib")
install.packages("dendextend")

library("DESeq2")
library('ggplot2') 
library('geneplotter')
library('vsn')
library('properties')

install.packages("rafalib")
install.packages("dendextend")

#get input values
args <- commandArgs(trailingOnly = TRUE)
algorithm <-args[1]
config_file <- args[2]
script_dir <-args[3]


function_script = paste(script_dir, "functions.R", sep='/')
source(function_script)

config <- read.properties(config_file)
dirlist = get_dir_info(algorithm, config, packagename="Deseq")
input.folder = dirlist$input.dir
output.folder= dirlist$out.dir
count.matrix.file = config$matrix_file
metadata.file = config$metadata_file
specifiedFactor = config$cluster_factor
 
#------------------------ 1. loading count matrix file ----------------
sprintf("--1. loading count matirx file %s",count.matrix.file)

if (endsWith(count.matrix.file, "gz")) {
  cts <- read.table(gzfile(count.matrix.file), sep="\t", header=TRUE)
} else {
  cts <- read.table(count.matrix.file, sep="\t", header=TRUE)
}
print(cts)

#----------------------- 2. filtering data - only keep 2000 genes based on variance ---------
print("-- 2. Trim Gene by Variance")
library(countSubsetNorm)
cts=var_filter_counts(cts, n_keep=2000, genes="rows",log2_transform=FALSE)

#-----------------------3. finding factors defined in metadata ---------------------------
sprintf("--3. reading all factors defined in metatdata %s",metadata.file )

factor_info=get_factors_info(metadata.file)
fctrs=factor_info$factors
factorlists=factor_info$factorlists

colData=as.data.frame(factorlists)
colnames(colData)=fctrs
rownames(colData)=colnames(cts)
print("-- feature data frame created")

designvar=paste(fctrs, collapse =" + ")
designvar=as.formula(paste('~',designvar,sep=""))

#--------------------- 4. Creating DESeqDataSetFromMatrix ------------------------
print("---4, creating DESeqDataSetFromMatrix")
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = colData,
                              design= designvar)
 
print("--DESeqDataSet created")

##  pre-filter low count genes
keep <- rowSums(counts(dds)) >= 25
dds <- dds[keep,]

#Extracting transformed values
vsd <- vst(dds, blind=FALSE)
print("--vst completed")
head(assay(vsd), 3)

#Ploting SD
sdPlot = paste(output.folder, '/meanSd.png', sep="")
png(sdPlot)
meanSdPlot(assay(vsd))
dev.off()

#---------------------- 5. Performing counts normalization based on the input algorithm ------------
print("--5. Performing normalization")
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
  
  mdsplot = paste(output.folder, '/mds.png', sep="")
  png(mdsplot)
  plotMDS.default(dgList) 
  dev.off()
  
  sizeFactors(dds)=dgList$samples$norm.factors
}  
print("--Normalization complete")

#---------------------6. Plotting --------------------------------------
print("-- Start to generate normalization plots")
normPlot1 = paste(output.folder, '/normalization.png', sep="")
png(normPlot1)
par(mfrow=c(1,1))
plot(log2( 1 + counts(dds, normalized=TRUE)[ , 1:2] ),col=rgb(0,0,0,.2), pch=16, cex=0.3 )
dev.off() 

normPlot2 = paste(output.folder, '/normalization_raw.png', sep="")
png(normPlot2)
plot( log2( 1 + counts(dds)[ , 1:2] ), col=rgb(0,0,0,.2), pch=16, cex=0.3 )
dev.off()

normPlot3 = paste(output.folder, '/normalization_vsd.png', sep="")
png(normPlot3)
plot( assay(vsd)[ , 1:2], col=rgb(0,0,0,.2), pch=16, cex=0.3 )
dev.off()
 

#--plotting meanSd
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

#--- plotting pheatmap
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=FALSE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,fctrs])
rownames(df)<-colnames(assay(vsd)[select,])
heatmap = paste(output.folder, '/pheatmap.png', sep="")
png(heatmap)
 
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
 
dev.off()
print("--pheatmap created")

#---Hierarchical cluster analysis

#-------------------plotting heatmap -------------
library(RColorBrewer) 
library("rafalib")
sampleDists <- dist(t(assay(vsd)))

sampleDistMatrix <- as.matrix( sampleDists )
coldf = colData(dds)
colindex <-grep(specifiedFactor, colnames(coldf))
rownames(sampleDistMatrix) <- paste(coldf[,colindex], sep="-" )

colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

cols=palette(brewer.pal(8, "Dark2"))[unlist(lapply(colData[specifiedFactor], as.numeric))]
hc <- hclust(sampleDists, method="ward.D")
 
library(gplots)
heatmap = paste(output.folder, '/heatmap_Eucl.png', sep="")
png(heatmap)
heatmap.2( sampleDistMatrix, Rowv=as.dendrogram(hc),
           symm=TRUE, trace="none", col=colors,
           ColSideColors = cols,
           margins=c(1,5), labCol=FALSE )
dev.off()

library(dendextend)
dendro=paste(output.folder, '/heatmap_dendro_Eucl.png', sep="")
dend <- as.dendrogram(hc)
dend1 <- color_branches(dend, k = 4)
colors_to_use <- unlist(lapply(colData[specifiedFactor], as.numeric))
colors_to_use <- colors_to_use[order.dendrogram(dend)]
labels_colors(dend1) <- colors_to_use
png(dendro)
plot(dend1, main = "Euclidean Distance Clustering")
dev.off()

#---------------- spearman heatmap ------------

library(psych)
sampleSpearman=cor(assay(vsd),method="spearman")
sampleSpearmanDist=cor2dist(sampleSpearman)

colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
plot(hclust(as.dist(sampleSpearmanDist), method="ward.D2"))
cols=palette(brewer.pal(8, "Dark2"))[unlist(lapply(colData[specifiedFactor], as.numeric))]
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
dendro1=paste(output.folder, '/heatmapSpearman_dendro.png', sep="")
dend <- as.dendrogram(hc)
dend1 <- color_branches(dend, k = 4)
colors_to_use <- unlist(lapply(colData[specifiedFactor], as.numeric))
colors_to_use <- colors_to_use[order.dendrogram(dend)]
labels_colors(dend1) <- colors_to_use
png(dendro1)
plot(dend1, main = "Spearman's rho Clustering")
dev.off()

#---------------pearson heatmap ------------------
samplePearson=cor(assay(vsd),method="pearson")
samplePearsonDist=cor2dist(samplePearson)
 
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)


cols=palette(brewer.pal(8, "Dark2"))[unlist(lapply(colData[specifiedFactor], as.numeric))]
hc <- hclust(as.dist(samplePearsonDist), method="ward.D")
#png(dendro_pearson)
#plot(hc)
#clusters=cutree(hc,k=4)

library(gplots)
heatmap = paste(output.folder, '/heatmap_pearson.png', sep="")
png(heatmap)
heatmap.2(samplePearsonDist, Rowv=as.dendrogram(hc),
          symm=TRUE, trace="none", col=colors,
          ColSideColors = cols,
          margins=c(1,5), labCol=FALSE )
dev.off()

print("--Heatmaps created")

#----------------Plotting PCA ------------------------
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

PCA1 = paste(output.folder, '/PCA_top_100.png', sep="")
png(PCA1, width=800, height=800)
plotPCA(vsd, intgroup=fctrs, ntop=100)
dev.off()

print("--PCA Completed")

dds <- DESeq(dds, test="LRT", reduced=~1, sfType="poscounts", useT=TRUE, minmu=1e-6, minReplicatesForReplace=Inf)
res <- results(dds)
resultsNames(dds)
print("DESeq completed")

#write summary to a file
summary.file = paste(output.folder, "resSum.txt" ,sep="/" )
sink(summary.file)
summary(res)
sink()
 

resOrdered <- res[order(res$pvalue),]
head(resOrdered)
output.file = paste(output.folder, "resOrdered.txt" ,sep="/" )
#Write DESeq results to file (no gene names)
write.table(resOrdered,output.file,row.names = TRUE, col.names = TRUE)

#--shrinks log-fold change values for plotMA
library("apeglm")
resLFC <- lfcShrink(dds, coef=resultsNames(dds)[2], type="apeglm")
resLFC

#--outlier detection using cooks distance
random=base::sample(x=ncol(assays(dds)[["cooks"]]),size=10,replace=FALSE)

cooks=paste(output.folder, '/cooks_distance.png', sep="")
png(cooks)
boxplot(log10(assays(dds)[["cooks"]])[,random], range=0, las=2)
dev.off()
#---------------------------------------------------------------------------------
#----------  cell-cell clustering based on expression of top 100 HVG's
#----------------------------------------------------------------------------------

#--------heatmap clustered by euclidean distance
sampleDists <- dist(t(assay(vsd[rownames(resOrdered)[1:100],])))
plot(hclust(sampleDists, method="ward.D"))
sampleDistMatrix <- as.matrix( sampleDists )

coldf = colData(dds)
colindex <-grep(specifiedFactor, colnames(coldf))
rownames(sampleDistMatrix) <- paste(coldf[,colindex], sep="-" )

colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
library("rafalib")
cols=palette(brewer.pal(8, "Dark2"))[unlist(lapply(colData[specifiedFactor], as.numeric))]
hc <- hclust(sampleDists, method="ward.D")
library(gplots)

heatmap = paste(output.folder, '/heatmap_Eucl_100.png', sep="")
png(heatmap)
heatmap.2( sampleDistMatrix, Rowv=as.dendrogram(hc),
           symm=TRUE, trace="none", col=colors,
           ColSideColors = cols,
           margins=c(1,5), labCol=FALSE )
dev.off()

library(dendextend)
dend <- as.dendrogram(hc)
dend1 <- color_branches(dend, k = 4)
colors_to_use <- unlist(lapply(colData[specifiedFactor], as.numeric))
colors_to_use <- colors_to_use[order.dendrogram(dend)]
labels_colors(dend1) <- colors_to_use
dendro2 = paste(output.folder, '/dendro_Eucl_100.png', sep="")
png(dendro2)
plot(dend1, main = "Euclidean Distance Clustering using top 100 HVGs")
dev.off()

#--------Spearman correlation heatmap ----------------------
library(psych)
sampleSpearman=cor(assay(vsd)[rownames(resOrdered)[1:100],],method="spearman")
sampleSpearmanDist=cor2dist(sampleSpearman)

rownames(sampleSpearmanDist) <- paste(vsd$molecular_subtype, sep="-" )
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
#plot(hclust(as.dist(sampleSpearmanDist),method="ward.D"))

cols=palette(brewer.pal(8, "Dark2"))[unlist(lapply(colData[specifiedFactor], as.numeric))]
hc <- hclust(as.dist(sampleSpearmanDist), method="ward.D")
library(gplots)
heatmap = paste(output.folder, '/heatmap_Spearman_100.png', sep="")
png(heatmap)
heatmap.2(sampleSpearmanDist, Rowv=as.dendrogram(hc),
          symm=TRUE, trace="none", col=colors,
          ColSideColors = cols,
          margins=c(1,5), labCol=FALSE )


library(dendextend)

dend <- as.dendrogram(hc)
dend1 <- color_branches(dend, k = 4)
colors_to_use <- unlist(lapply(colData[specifiedFactor], as.numeric))
colors_to_use <- colors_to_use[order.dendrogram(dend)]
labels_colors(dend1) <- colors_to_use
dendro2 = paste(output.folder, '/dendro_Spearman_100.png', sep="")
png(dendro2)
plot(dend1, main = "Spearman rho Clustering using top 100 HVGs")
dev.off()
#--------Pearson correlation heatmap ----------------------
samplePearson=cor(assay(vsd)[rownames(resOrdered)[1:100],],method="pearson")
samplePearsonDist=cor2dist(samplePearson)

rownames(samplePearsonDist) <- paste(vsd$molecular_subtype, sep="-" )
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

cols=palette(brewer.pal(8, "Dark2"))[unlist(lapply(colData[specifiedFactor], as.numeric))]
hc <- hclust(as.dist(samplePearsonDist), method="ward.D")

#plot(hc)
#clusters=cutree(hc,k=4)

library(gplots)
heatmap = paste(output.folder, '/heatmap_Pearson_100.png', sep="")
png(heatmap)
heatmap.2(samplePearsonDist, Rowv=as.dendrogram(hc),
          symm=TRUE, trace="none", col=colors,
          ColSideColors = cols,
          margins=c(1,5), labCol=FALSE )
dev.off()
 
library(dendextend)

dend <- as.dendrogram(hc)  
dend1 <- color_branches(dend, k = 4)
colors_to_use <- unlist(lapply(colData[specifiedFactor], as.numeric))
colors_to_use <- colors_to_use[order.dendrogram(dend)]
labels_colors(dend1) <- colors_to_use
dendro2 = paste(output.folder, '/dendro_Pearson_100.png', sep="")
png(dendro2)
plot(dend1, main = "Pearson rho Clustering with top 100 HVGs")
dev.off()
print("--Heatmaps created")

#---------Plotting PCA ----------

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


MA = paste(output.folder, '/MA.png', sep="")
png(MA)
plotMA(resLFC, ylim=c(-10,20))
dev.off()
#------Top Genes
topGene <- rownames(res)[which.min(res$padj)]

par( mfrow = c( 1, 1 ) )
cts = paste(output.folder, '/counts.png', sep="")
png(cts)
plotCounts(dds, gene=topGene, intgroup=fctrs)
dev.off()

MA1 = paste(output.folder, '/MA1.png', sep="")
png(MA1)
DESeq2::plotMA(res, ylim=c(-10,20))
with(res[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})
dev.off()

DispEsts = paste(output.folder, '/DispEsts.png', sep="")
png(DispEsts)
plotDispEsts(dds)
dev.off()

hist1 = paste(output.folder, '/hist1.png', sep="")
png(hist1)
hist(res$pvalue, breaks=20, col="grey50", border="white")
dev.off()

hist2 = paste(output.folder, '/hist2.png', sep="")
png(hist2)
hist(res$pvalue[res$baseMean > 1], breaks=20, col="grey50", border="white")
dev.off()

#----- Takes most highly variable genes
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

#--add gene names to result matrix
library("AnnotationDbi")
library("org.Hs.eg.db")
library("org.Mm.eg.db") 

tmp=gsub("\\..*","",row.names(res)) 
 
genome=config$genome_type
if ( genome=="human" ) {
  genomeAnnotation=org.Hs.eg.db
}else if ( genome=="mouse" ){
  genomeAnnotation=org.Mm.eg.db
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

output.file = paste(output.folder, "resOrdered_final.txt" ,sep="/" )
write.table(resOrdered,output.file,row.names = TRUE, col.names = TRUE)

