library(DESeq2)
library("ggplot2") 
library("geneplotter")
library('vsn') 

#pass in arguments
args <- commandArgs(trailingOnly = TRUE) 

algorithm <- args[1]
input.folder <- args[2]
print(input.folder) 
output.folder <- args[3]
genomeAnnotation <- args[4]

metadata.file<- args[5]

specifiedFactor<-args[6]

#convert to lowercase
algorithm <- tolower(algorithm) 
algorithm <- "default"


# read input files
setwd(input.folder) 
filelist = list.files(pattern="*counts.txt")
print(filelist)

#assuming tab separated values with a header   
datalist = lapply(filelist, function(x)read.table( x, header=TRUE, row.names=1, fill = FALSE)) 
cts = do.call("cbind", datalist) 
#remove last 5 rows which contains summary info

cts <- cts[1:(nrow(cts)-5),]



#trim genes based on variance
library(countSubsetNorm)
cts=var_filter_counts(cts, n_keep=10000, genes="rows",log2_transform=FALSE)
print("--Genes trimmed by Variance")

#read metadata.file
df <- read.table(metadata.file, header = TRUE, sep = "\t")
fctrs=colnames(df)

#create colData matrix
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
#colData=df
colData=as.data.frame(factorlists)
colnames(colData)=fctrs
rownames(colData)=colnames(cts)

print("colData created")

designvar=paste(fctrs, collapse =" + ")
designvar=as.formula(paste('~',designvar,sep=""))

#make deseq data set
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = colData,
                              design= designvar)

print("--DESeqDataSet created")
keep <- rowSums(counts(dds)) >= 25
dds <- dds[keep,]

#transform dataset counts for visuals
vsd <- vst(dds, blind=FALSE)
print("--vst completed")
head(assay(vsd), 3)


#sdPlot = paste(output.folder, '/meanSd.png', sep="")
#png(sdPlot)
#meanSdPlot(assay(vsd))
#dev.off()

#normalization of counts
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
  gene_len <-read.delim(paste(input.folder,"hg38_gene_length.txt",sep=""), header=F, sep="\t")
  colnames(gene_len) <- c("GeneName","Len")   
  mcols(dds)$basepairs <- gene_len$Len
  fpkm(dds, gene.lengths)  
}

#normalized counts plot
print("--Normalization complete")
normPlot1 = paste(output.folder, '/normalization.png', sep="")
png(normPlot1)
par(mfrow=c(1,1))
plot( log2( 1 + counts(dds, normalized=TRUE)[ , 1:2] ),col=rgb(0,0,0,.2), pch=16, cex=0.3 )
dev.off()

#raw counts plot
normPlot2 = paste(output.folder, '/normalization_raw.png', sep="")
png(normPlot2)
plot( log2( 1 + counts(dds)[ , 1:2] ), col=rgb(0,0,0,.2), pch=16, cex=0.3 )
dev.off()

#vsd plot
normPlot3 = paste(output.folder, '/normalization_vsd.png', sep="")
png(normPlot3)
plot( assay(vsd)[ , 1:2], col=rgb(0,0,0,.2), pch=16, cex=0.3 ) 

dev.off()

notAllZero <- (rowSums(counts(dds))>0)

# Non-Zero
random=base::sample(x=ncol(counts(dds)),size=25,replace=FALSE)

#normalized counts meansd plot
meanSd1 = paste(output.folder, '/meanSd_norm.png', sep="")
png(meanSd1)
meanSdPlot(log2(counts(dds,normalized=TRUE)[notAllZero,] + 1)) 
dev.off()

#raw counts meansd plot
meanSd2 = paste(output.folder, '/meanSd_raw.png', sep="")
png(meanSd2)
meanSdPlot(log2(counts(dds)[notAllZero,] + 1))
dev.off()

#vsd meansd plot
meanSd3 = paste(output.folder, '/meanSd_vsd.png', sep="")
png(meanSd3)
meanSdPlot(assay(vsd[notAllZero,]))
dev.off()

#box plot of raw counts from a random number of samples
boxplot = paste(output.folder, '/boxplot_raw.png', sep="")
png(boxplot)
boxplot(log2(counts(dds)[notAllZero,random] + 1)) 
dev.off()

#box plot of normalized counts from a random number of samples
boxplot2 = paste(output.folder, '/boxplot_norm.png', sep="")
png(boxplot2)
boxplot(log2(counts(dds,normalized=TRUE)[notAllZero,random] + 1)) 
dev.off()

#heatmap comparison between cells and highly expressed genes
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=FALSE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds))
rownames(df)<-colnames(assay(vsd)[select,])
heatmapp = paste(output.folder, '/pheatmap.png', sep="")
png(heatmapp)
pheatmap(assay(vsd)[select,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=df, clustering_method = "ward.D", show_colnames = FALSE)
dev.off()



print("--pheatmap created")

#cell to cell heatmap (Euclidean distance)

sampleDists <- dist(t(assay(vsd)))
plot(hclust(sampleDists,method="ward.D")) 

library(RColorBrewer)

sampleDistMatrix <- as.matrix( sampleDists )

rownames(sampleDistMatrix) <- paste(colData[,specifiedFactor], sep="-" )
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
library("rafalib")
cols=palette(brewer.pal(8, "Dark2"))[as.fumeric(as.character(colData[,specifiedFactor]))]
hc <- hclust(sampleDists, method="ward.D")
library(gplots)
heatmap = paste(output.folder, '/heatmapEuclidean.png', sep="")
png(heatmap)
heatmap.2( sampleDistMatrix, Rowv=as.dendrogram(hc),
           symm=TRUE, trace="none", col=colors,
           ColSideColors = cols,
           margins=c(1,5), labCol=FALSE ,
           main = "Cell-to-cell Euclidean distance")
dev.off()

#Euclidean distance dendrogram
library(dendextend)

dend <- as.dendrogram(hc)
# Like: 
# dend <- USArrests[1:5,] %>% dist %>% hclust %>% as.dendrogram

dend1 <- color_branches(dend, k = 3)
colors_to_use <- as.numeric(colData[,specifiedFactor])
colors_to_use <- colors_to_use[order.dendrogram(dend)]
labels_colors(dend1) <- colors_to_use

dendr = paste(output.folder, '/dendrEuclidean.png', sep="")
png(dendr)
plot(dend1, main = "Euclidean Distance Clustering")
dev.off()

#cell to cell heatmap (Spearman correlation)

library(psych)
sampleSpearman=cor(assay(vsd),method="spearman")
sampleSpearmanDist=cor2dist(sampleSpearman)

rownames(sampleSpearmanDist) <- paste(colData[,specifiedFactor], sep="-" )
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
plot(hclust(as.dist(sampleSpearmanDist), method="ward.D2"))
cols=palette(brewer.pal(8, "Dark2"))[as.fumeric(as.character(colData[,specifiedFactor]))]
hc <- hclust(as.dist(sampleSpearmanDist), method="ward.D")
library(gplots)
heatmapS = paste(output.folder, '/heatmapSpearman.png', sep="")
png(heatmapS)
heatmap.2(sampleSpearmanDist, Rowv=as.dendrogram(hc),
          symm=TRUE, trace="none", col=colors,
          ColSideColors = cols,
          margins=c(1,5), labCol=FALSE , labRow = FALSE ,
          main = "Cell-to-cell Spearman distance")
dev.off()

#spearman distance dendrogram
library(dendextend)

dend <- as.dendrogram(hc)
# Like: 
# dend <- USArrests[1:5,] %>% dist %>% hclust %>% as.dendrogram

dend1 <- color_branches(dend, k = 4)
colors_to_use <- as.numeric(colData[,specifiedFactor])
colors_to_use <- colors_to_use[order.dendrogram(dend)]
labels_colors(dend1) <- colors_to_use



dendr1 = paste(output.folder, '/dendrSpearman.png', sep="")
png(dendr1)
plot(dend1, main = "Spearman Distance Clustering")
dev.off()

#cell-to-cell heatmap (Pearson correlation)
samplePearson=cor(assay(vsd),method="pearson")
samplePearsonDist=cor2dist(samplePearson)
plot(hclust(as.dist(samplePearsonDist),method="ward.D"))
rownames(samplePearsonDist) <- paste(colData[,specifiedFactor], sep="-" )
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

cols=palette(brewer.pal(8, "Dark2"))[as.fumeric(as.character(colData[,specifiedFactor]))]
hc <- hclust(as.dist(samplePearsonDist), method="ward.D")
plot(hc)
#clusters=cutree(hc,k=4)

library(gplots)
heatmapP = paste(output.folder, '/heatmapPearson.png', sep="")
png(heatmapP)
heatmap.2(samplePearsonDist, Rowv=as.dendrogram(hc),
          symm=TRUE, trace="none", col=colors,
          ColSideColors = cols,
          margins=c(1,5), labCol=FALSE , labRow = FALSE ,
          main = "Cell-to-cell Pearson distance")
dev.off()

print("--Heatmaps created")

#Pearson distance dendrogram
library(dendextend)

dend <- as.dendrogram(hc)
# Like: 
# dend <- USArrests[1:5,] %>% dist %>% hclust %>% as.dendrogram

dend1 <- color_branches(dend, k = 4)
colors_to_use <- as.numeric(colData[,specifiedFactor])
colors_to_use <- colors_to_use[order.dendrogram(dend)]
labels_colors(dend1) <- colors_to_use

dendr = paste(output.folder, '/dendrPearson.png', sep="")
png(dendr)
plot(dend1, main = "Pearson Distance Clustering")
dev.off()

#PCA
PCA = paste(output.folder, '/PCA.png', sep="")
png(PCA, width=800, height=800)
par(mfrow=c(1,1))
plotPCA(vsd, intgroup=specifiedFactor)+ coord_fixed(ratio=3) 

dev.off()

print("--PCA Completed")

#DESeq analysis
dds <- DESeq(dds, test="LRT", reduced=~1, sfType="poscounts", useT=TRUE, minmu=1e-6, minReplicatesForReplace=Inf)
res <- results(dds)

print("DESeq completed")

#library("IHW")
#resIHW <- results(dds, filterFun=ihw)
#summary(resIHW)
 
resultsNames(dds)

#lfcshrink
resLFC <- lfcShrink(dds, coef=resultsNames(dds)[2], type="apeglm")
resLFC

summary(res)

#sort results by pvalue and annotate
library("AnnotationDbi")


#convertIDs <- function( ids, from, to, db, ifMultiple=c("putNA", "useFirst")) {
#  stopifnot( inherits( db, "AnnotationDb" ) )
#  ifMultiple <- match.arg( ifMultiple )
#  suppressWarnings( selRes <- AnnotationDbi::select(
#    db, keys=ids, keytype=from, columns=c(from,to) ) )
#  if ( ifMultiple == "putNA" ) {
#    duplicatedIds <- selRes[ duplicated( selRes[,1] ), 1 ]
#    selRes <- selRes[ ! selRes[,1] %in% duplicatedIds, ]
#  }
#  return( selRes[ match( ids, selRes[,1] ), 2 ] )
#}

tmp=gsub("\\..*","",row.names(res))
if ( genome=="human" ) {
  library("org.Hs.eg.db")
  genomeAnnotation= org.Hs.eg.db
}else if ( genome=="mouse" ){
  library("org.Mm.eg.db")
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

#sort results by p value
resOrdered <- res[order(res$pvalue),]
head(resOrdered)
resOrdered=na.omit(resOrdered)
resOrdered=resOrdered[rownames(resOrdered) %in% rownames(assay(vsd)),]
write.table(resOrdered,"resOrdered.txt",row.names = TRUE, col.names = TRUE)
head(resOrdered)

#outlier check
#random=base::sample(x=ncol(assays(dds)[["cooks"]]),size=10,replace=FALSE)
#boxplot(log10(assays(dds)[["cooks"]])[,random], range=0, las=2)



#library("Rfast")
#cvmatrix=rowcvs(counts(dds))
#cvgenematrix=cbind(rownames(counts(dds)),cvmatrix)
#cvOrdered=cvgenematrix[order(cvgenematrix[,2],decreasing = TRUE)]

#cell-to-cell heatmap, filtered genes (Euclidean Distance)
sampleDists <- dist(t(assay(vsd[rownames(resOrdered)[1:1000],])))
plot(hclust(sampleDists, method="ward.D")) 
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(colData[,specifiedFactor], sep="-" )
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
library("rafalib")
cols=palette(brewer.pal(8, "Dark2"))[as.fumeric(as.character(colData[,specifiedFactor]))]
hc <- hclust(sampleDists, method="ward.D")
library(gplots)
heatmap = paste(output.folder, '/heatmapRes.png', sep="")
png(heatmap)
heatmap.2( sampleDistMatrix, Rowv=as.dendrogram(hc),
           symm=TRUE, trace="none", col=colors,
           ColSideColors = cols,
           margins=c(1,5), labCol=FALSE, labRow = FALSE ,
           main = "Cell-to-cell Euclidean distance (lowest 1000 genes by p-value)" )
dev.off()

#Euclidean distance dendrogram
library(dendextend)

dend <- as.dendrogram(hc)

dend1 <- color_branches(dend, k = 4)
colors_to_use <- as.numeric(colData[,specifiedFactor])
colors_to_use <- colors_to_use[order.dendrogram(dend)]
labels_colors(dend1) <- colors_to_use
dendr = paste(output.folder, '/dendrResEuclidean.png', sep="")
png(dendr)
plot(dend1, main = "Euclidean Distance Clustering, genes trimmed by p-value (lowest 1000)")
dev.off()

#cell-to-cell heatmap, filtered genes (Spearman Correlation)
library(psych)
sampleSpearman=cor(assay(vsd)[rownames(resOrdered)[1:1000],],method="spearman")
sampleSpearmanDist=cor2dist(sampleSpearman)

rownames(sampleSpearmanDist) <- paste(colData[,specifiedFactor], sep="-" )
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
plot(hclust(as.dist(sampleSpearmanDist),method="ward.D"))
cols=palette(brewer.pal(8, "Dark2"))[as.fumeric(as.character(colData[,specifiedFactor]))]
hc <- hclust(as.dist(sampleSpearmanDist), method="ward.D")
library(gplots)
heatmap = paste(output.folder, '/heatmapResS.png', sep="")
png(heatmap)
heatmap.2(sampleSpearmanDist, Rowv=as.dendrogram(hc),
          symm=TRUE, trace="none", col=colors,
          ColSideColors = cols,
          margins=c(1,5), labCol=FALSE ,labRow = FALSE ,
          main = "Cell-to-cell Spearman distance (lowest 1000 genes by p-value)")

library(dendextend)


#Spearman distance dendrogram


dend1 <- color_branches(dend, k = 4)
colors_to_use <- as.numeric(colData[,specifiedFactor])
colors_to_use <- colors_to_use[order.dendrogram(dend)]
labels_colors(dend1) <- colors_to_use

dendr = paste(output.folder, '/dendrResSpearman.png', sep="")
png(dendr)
plot(dend1, main = "Spearman Distance Clustering, genes trimmed by p-value (lowest 1000)")
dev.off()



#cell-to-cell heatmap, filtered genes (Pearson Correlation)
samplePearson=cor(assay(vsd)[rownames(resOrdered)[1:1000],],method="pearson")
samplePearsonDist=cor2dist(samplePearson)
plot(hclust(as.dist(samplePearsonDist),method="ward.D"))

rownames(samplePearsonDist) <- paste(colData[,specifiedFactor], sep="-" )
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

cols=palette(brewer.pal(8, "Dark2"))[as.fumeric(as.character(colData[,specifiedFactor]))]
hc <- hclust(as.dist(samplePearsonDist), method="ward.D")
plot(hc)

library(gplots)
heatmap = paste(output.folder, '/heatmapResP.png', sep="")
png(heatmap)
heatmap.2(samplePearsonDist, Rowv=as.dendrogram(hc),
          symm=TRUE, trace="none", col=colors,
          ColSideColors = cols,
          margins=c(1,5), labCol=FALSE, labRow = FALSE ,
          main = "Cell-to-cell Pearson distance (lowest 1000 genes by p-value)")
dev.off()


# install.packages("dendextend")
library(dendextend)

dend <- as.dendrogram(hc)
# Like: 
# dend <- USArrests[1:5,] %>% dist %>% hclust %>% as.dendrogram

dend1 <- color_branches(dend, k = 4)
colors_to_use <- as.numeric(colData[,specifiedFactor])
colors_to_use <- colors_to_use[order.dendrogram(dend)]
labels_colors(dend1) <- colors_to_use

#Pearson distance dendrogram
dendr = paste(output.folder, '/dendrResPearson.png', sep="")
png(dendr)
plot(dend1, main = "Pearson Distance Clustering, genes trimmed by p-value (lowest 1000)")
dev.off()



print("--Heatmaps created")

#PCA, genes filtered
PCAres = paste(output.folder, '/PCAres.png', sep="")
png(PCAres, width=800, height=800)
par(mfrow=c(1,1))
plotPCA(vsd[rownames(resOrdered)[1:1000]], 
        intgroup=specifiedFactor)+coord_fixed(ratio=3)+ ggtitle("PCA for genes filtered by p-value (lowest 1000)")

dev.off()

#MA plot
plotDispEsts(dds)

plotMA(resLFC, ylim=c(-20,20))
       
#cell-to-gene heatmap
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
sidecols <- palette(brewer.pal(8, "Dark2"))[ colData[,specifiedFactor] ]
mat <- assay(vsd)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
heatmapTVG = paste(output.folder, '/heatmapTVG.png', sep="")
png(heatmapTVG)
heatmap.2(mat, trace="none", col=colors, ColSideColors=sidecols,
          labRow=FALSE, mar=c(10,2), scale="row", Colv=as.dendrogram(hc))
dev.off()

#TSNE
library(ggplot2)
library(Rtsne)
tsne=Rtsne(t(assay(vsd)), perplexity=5)

tsneplot = paste(output.folder, '/tsne5.png', sep="")
png(tsneplot)
tsne_plot <- data.frame(x = tsne$Y[,1], y = tsne$Y[,2], col = colData[,specifiedFactor])
ggplot(tsne_plot) + geom_point(aes(x=x, y=y, color=col))+ggtitle("t-SNE for unfiltered genes, Perplexity = 5")
dev.off()

#TSNE, genes filtered
tsneRes=Rtsne(t(assay(vsd[rownames(resOrdered)[1:1000]])), perplexity=5)

tsneResPlot = paste(output.folder, '/tsne5resOrdered.png', sep="")
png(tsneResPlot)
tsne_plot <- data.frame(x = tsneRes$Y[,1], y = tsneRes$Y[,2], col = colData[,specifiedFactor])
ggplot(tsne_plot) + geom_point(aes(x=x, y=y, color=col))+ggtitle("t-SNE for genes filtered by p-value (lowest 1000), Perplexity = 5")
dev.off()

library(Seurat)

for ( i in levels(colData(dds)[,specifiedFactor]) ) {
  markers <- FindMarkers(counts(dds), cells.1=rownames(colData(dds)[colData(dds)[,specifiedFactor]==i,]), 
    cells.2=rownames(colData(dds)[!colData(dds)[,specifiedFactor]==i,]), test.use = "DESeq2", slot = "counts")
  write.table(markers,file=paste(gsub(" ", "", as.character(i), fixed = TRUE),".txt", sep=""),".txt", col.names=1)

}

markers <- FindMarkers(counts(dds), cells.1=rownames(colData(dds)[colData(dds)$molecular_subtype=="triple-negative breast cancer (TNBC)",]), 
                       cells.2=rownames(colData(dds)[!colData(dds)$molecular_subtype=="triple-negative breast cancer (TNBC)",]), test.use = "DESeq2", slot = "counts")
head(x = markers)

markers1 <- FindMarkers(counts(dds), cells.1=rownames(colData(dds)[colData(dds)$molecular_subtype=="estrogen receptor positive (ER+)",]), 
                       cells.2=rownames(colData(dds)[!colData(dds)$molecular_subtype=="estrogen receptor positive (ER+)",]), test.use = "DESeq2",slot = "counts")

markers2 <- FindMarkers(counts(dds), cells.1=rownames(colData(dds)[colData(dds)$molecular_subtype=="double positive (ER+ and HER2+)",]), 
                        cells.2=rownames(colData(dds)[!colData(dds)$molecular_subtype=="double positive (ER+ and HER2+)",]), test.use = "DESeq2",slot = "counts")

markers3 <- FindMarkers(counts(dds), cells.1=rownames(colData(dds)[colData(dds)$molecular_subtype=="human epidermal growth factor receptor 2 positive (HER2+)",]), 
                        cells.2=rownames(colData(dds)[!colData(dds)$molecular_subtype=="human epidermal growth factor receptor 2 positive (HER2+)",]), test.use = "DESeq2",slot = "counts")
