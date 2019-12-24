library(R.utils)
library(gdata)
library(scater)

args <- commandArgs(trailingOnly = TRUE)

algorithm<- args[1]
input.folder<- args[2]
output.folder<- args[3]
genomeAnnotation<-args[4]
metadata.file<-args[5]
#convert to lowercase
algorithm <- tolower(algorithm)

readFormat <- function(infile) {
  # First column is empty.
  metadata <- read.delim(infile, stringsAsFactors=FALSE, header=FALSE, nrow=10)[,-1]
  rownames(metadata) <- metadata[,1]
  metadata <- metadata[,-1]
  metadata <- as.data.frame(t(metadata))
  # First column after row names is some useless filler.
  counts <- read.delim(infile, stringsAsFactors=FALSE, header=FALSE, row.names=1, skip=11)[,-1]
  counts <- as.matrix(counts)
  return(list(metadata=metadata, counts=counts))
}

#endo.data <- readFormat("test.txt")

# read input files
setwd(input.folder)
filelist = list.files(pattern="*counts.txt")
gene_len <-read.delim(paste(input.folder,"hg38_gene_length.txt",sep=""), header=F, sep="\t")
colnames(gene_len) <- c("GeneName","Len")

datalist = lapply(filelist, function(x)read.table( x, header=TRUE, row.names=1, fill = FALSE))
cts = do.call("cbind", datalist)

#cts=read.table("./counts.txt")

#remove last 5 rows which contains summary info
cts <-cts[1:(nrow(cts)-5),]

df <- read.table(metadata.file, header = TRUE, sep = "\t")
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
print("colData created")


sce <- SingleCellExperiment(assays=list(counts=as.matrix(cts)), colData=as.matrix(colData))
sce$total_features
dim(sce)
print("--SingleCellExperiment Created")
 
is.spike <- grepl("^ERCC", rownames(sce))
is.mito <- grepl("^mt-", rownames(sce))

sce <- calculateQCMetrics(sce, feature_controls=list(ERCC=is.spike, Mt=is.mito))
#head(colnames(pData(sce)))

library(scran)
#isSpike(sce) <- "ERCC"

par(mfrow=c(1,2))

#hist(sce$total_counts/1e6, xlab="Library sizes (millions)", main="", breaks = 20, col="grey80", ylab="Number of cells")
#hist(as.numeric(sce$total_features_by_counts), xlab="Number of expressed genes", main="", breaks=20, col="grey80", ylab="Number of cells")


#libsize.drop <- isOutlier(sce$total_counts, nmads=3, type="lower", log=TRUE)
#feature.drop <- isOutlier(sce$total_features_by_counts, nmads=3, type="lower", log=TRUE)
histogram = paste(output.folder, '/histogram.png', sep="")
png(histogram)
par(mfrow=c(1,2))
hist(sce$total_counts/1e6, xlab="Library sizes (millions)", main="",breaks=20, col="grey80", ylab="Number of cells")
hist(sce$total_features_by_counts, xlab="Number of expressed genes", main="",breaks=20, col="grey80", ylab="Number of cells")
dev.off()

libsize.drop <- isOutlier(sce$total_counts, nmads=3, type="lower", log=TRUE)
feature.drop <- isOutlier(sce$total_features_by_counts, nmads=3, type="lower", log=TRUE)
print("--Cells filtered")


par(mfrow=c(1,2))
hist(as.numeric(as.character(sce$pct_counts_Mt)), xlab="Mitochondrial proportion (%)",
ylab="Number of cells", breaks=20, main="", col="grey80")
hist(sce$pct_counts_ERCC, xlab="ERCC proportion (%)",
ylab="Number of cells", breaks=20, main="", col="grey80")

mito.drop <- isOutlier(sce$pct_counts_Mt, nmads=3, type="higher")
spike.drop <- isOutlier(sce$pct_counts_ERCC, nmads=3, type="higher")

sce <- sce[,!(libsize.drop | feature.drop | mito.drop | spike.drop)]
data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop),
           ByMito=sum(mito.drop), BySpike=sum(spike.drop), Remaining=ncol(sce))

fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))
#plotPCA(sce, pca_data_input="pdata") + fontsize

mm.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))
#library(yaml)
#library(ensembldb)
#library(EnsDb.Hsapiens.v79)
#library(EnsDb.Mmusculus.v79)

#if (genome="human") {
#  genomeAnnotation= EnsDb.Hsapiens.v79
#}else if (genome="mouse"){
#  genomeAnnotation= EnsDb.Mmusculus.v79
#}


#tmp=gsub("\\..*","",row.names(sce))


#anno <- select(EnsDb.Hsapiens.v79,tmp,keytype='GENEID', column="SYMBOL")
#acensembl <- anno$GENEID[match(tmp, anno$GENEID)]
#assignments <- cyclone(sce, mm.pairs, gene.names=acensembl)
#cellcycle = paste(output.folder, '/cellcycle.png', sep="")
#png(cellcycle)
#plot(assignments$score$G1, assignments$score$G2M, xlab="G1 score", ylab="G2/M score", pch=16)
#dev.off()
#print(sce)
#print(assignments)
#print(assignments$phases)
#assignments$phases[is.na(assignments$phases)] <- "NA"
#sce <- sce[,assignments$phases=="G1"]
#print("--Cell Cycle Assignment complete")
ave.counts <- rowMeans(counts(sce))
keep <- ave.counts >= 1
sum(keep)

par(mfrow=c(1,1))

histogram1 = paste(output.folder, '/histogram1.png', sep="")
png(histogram1)
hist(log10(ave.counts), breaks=100, main="", col="grey80",
xlab=expression(Log[10]~"average count"))
abline(v=log10(1), col="blue", lwd=2, lty=2)

HiExprs = paste(output.folder, '/hiexprs.png', sep="")
png(HiExprs)
plotHighestExprs(sce, n=50) + fontsize
dev.off()

numcells <- nexprs(sce, byrow=TRUE)
alt.keep <- numcells >= 10
sum(alt.keep)

SmoothScatter = paste(output.folder, '/smoothScatter.png', sep="")
png(SmoothScatter)
smoothScatter(log10(ave.counts), numcells, xlab=expression(Log[10]~"average count"),
ylab="Number of expressing cells")
is.ercc <- isSpike(sce, type="ERCC")
points(log10(ave.counts[is.ercc]), numcells[is.ercc], col="red", pch=16, cex=0.5)
dev.off()

sce <- sce[keep,]

if (algorithm == "default") {
#computeSumFactors (deconvolution)

sce <- computeSumFactors(sce, sizes=c(5, 10, 15, 20))
summary(sizeFactors(sce))

SizeFactors = paste(output.folder, '/sizefactors.png', sep="")
png(SizeFactors)
plot(sizeFactors(sce), sce$total_counts/1e6, log="xy",
ylab="Library size (millions)", xlab="Size factor")
dev.off()
#sce <- computeSpikeFactors(sce, type="ERCC", general.use=FALSE)
sce <- normalize(sce)
} else if (algorithm == "cpm") {
# calculateCPM
endog_genes <- !rowData(sce)$is_feature_control
assay(sce, "normcounts") <- log2(calculateCPM(sce, use.size.factors = FALSE) + 1)
sce<- normalize(sce)

} else if (algorithm == "tmm") {
#NormalizeExprs(method="TMM")
endog_genes <- !rowData(sce)$is_feature_control
sce <- normaliseExprs(
  sce,
  method = "TMM",
  feature_set = endog_genes
)

} else if (algorithm == "fpkm") {
#calculateFPKM,
sce <- getBMFeatureAnnos(
  sce,
  filters = "ensembl_gene_id", 
  attributes = c(
    "ensembl_gene_id",
    "hgnc_symbol",
    "chromosome_name",
    "start_position",
    "end_position"
  ), 
  feature_symbol = "hgnc_symbol",
  feature_id = "ensembl_gene_id",
  biomart = "ENSEMBL_MART_ENSEMBL", 
  dataset = "hsapiens_gene_ensembl",
  host = "www.ensembl.org"
)
print("--Normalization Complete")

sce.ann <- sce[!is.na(rowData(sce)$ensembl_gene_id), ]
  eff_length <- 
  abs(rowData(sce.ann)$end_position - rowData(sce.ann)$start_position) / 1000
#plot(eff_length, rowMeans(counts(sce.ann)))
tpm(sce.ann) <- log2(calculateFPKM(sce.ann, eff_length) + 1)
}

ExpVar = paste(output.folder, '/expvar.png', sep="")
png(ExpVar)
plotExplanatoryVariables(sce, variables=c("molecular_subtype", "log10_total_features_by_counts","log10_total_counts","pct_counts_in_top_50_features")) + fontsize
dev.off()

var.fit <- trendVar(sce, use.spikes=FALSE)
var.out <- decomposeVar(sce, var.fit)

varout = paste(output.folder, '/varout.png', sep="")
png(varout)
plot(var.out$mean, var.out$total, pch=16, cex=0.6, xlab="Mean log-expression",
ylab="Variance of log-expression")
o <- order(var.out$mean)
lines(var.out$mean[o], var.out$tech[o], col="dodgerblue", lwd=2)
cur.spike <- isSpike(sce)
points(var.out$mean[cur.spike], var.out$total[cur.spike], col="red", pch=16)
dev.off()

hvg.out <- var.out[which(var.out$FDR <= 0.05 & var.out$bio >= 0.5),]
hvg.out <- hvg.out[order(hvg.out$bio, decreasing=TRUE),]
nrow(hvg.out)
write.table(file=paste(output.folder,"hsc_hvg.tsv", sep="/"), hvg.out, sep="\t", quote=FALSE, col.names=NA)
head(hvg.out)

plotexprs = paste(output.folder, '/plotexprs.png', sep="")
png(plotexprs)
plotExpression(sce, rownames(hvg.out)[1:10], jitter="jitter") + fontsize
dev.off()

set.seed(100)
var.cor <- correlatePairs(sce, subset.row=rownames(hvg.out))
write.table(file="hsc_cor.tsv", var.cor, sep="\t", quote=FALSE, row.names=FALSE)
head(var.cor)

sig.cor <- var.cor$FDR <= 0.05
summary(sig.cor)

library(limma)
adj.exprs <- exprs(sce)
adj.exprs <- removeBatchEffect(adj.exprs, batch=sce$molecular_subtype)
norm_exprs(sce) <- adj.exprs

library(RBGL)
g <- ftM2graphNEL(cbind(var.cor$gene1, var.cor$gene2)[sig.cor,],
                  W=NULL, V=NULL, edgemode="undirected")
cl <- highlyConnSG(g)$clusters
cl <- cl[order(lengths(cl), decreasing=TRUE)]
head(cl)

chosen <- unique(c(var.cor$gene1[sig.cor], var.cor$gene2[sig.cor]))
norm.exprs <- exprs(sce)[chosen,,drop=FALSE]
heat.vals <- norm.exprs - rowMeans(norm.exprs)
library(gplots)

chosen <- unique(c(var.cor$gene1[sig.cor], var.cor$gene2[sig.cor]))
top.hvg <- rownames(hvg.out)[1]

set.seed(100)

#heatmap(euclidean)
sampleDists <- dist(t(logcounts(sce)))
plot(hclust(sampleDists,method="ward.D")) 
library(RColorBrewer)
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(sce$molecular_subtype, sep="-" )
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
library("rafalib")
cols=palette(brewer.pal(8, "Dark2"))[as.data.frame(sce$molecular_subtype)[,1]]
hc <- hclust(sampleDists, method="ward.D")
library(gplots)
heatmap = paste(output.folder, '/heatmap.png', sep="")
png(heatmap)
heatmap.2( sampleDistMatrix, Rowv=as.dendrogram(hc),
           symm=TRUE, trace="none", col=colors,
           ColSideColors = cols,
           margins=c(1,5), labCol=FALSE )
dev.off()

#heatmap(spearman)
library(psych)
sampleSpearman=cor(logcounts(sce),method="spearman")
sampleSpearmanDist=cor2dist(sampleSpearman)

rownames(sampleSpearmanDist) <- paste(sce$molecular_subtype, sep="-" )
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
plot(hclust(as.dist(sampleSpearmanDist), method="ward.D2"))
cols=palette(brewer.pal(8, "Dark2"))[as.data.frame(sce$molecular_subtype)[,1]]
hc <- hclust(as.dist(sampleSpearmanDist), method="ward.D")
library(gplots)
heatmap = paste(output.folder, '/heatmapSpearman.png', sep="")
png(heatmap)
heatmap.2(sampleSpearmanDist, Rowv=as.dendrogram(hc),
          symm=TRUE, trace="none", col=colors,
          ColSideColors = cols,
          margins=c(1,5), labCol=FALSE )
dev.off()

samplePearson=cor(logcounts(sce),method="pearson")
samplePearsonDist=cor2dist(samplePearson)
plot(hclust(as.dist(samplePearsonDist),method="ward.D"))
rownames(samplePearsonDist) <- paste(sce$molecular_subtype, sep="-" )
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

cols=palette(brewer.pal(8, "Dark2"))[as.data.frame(sce$molecular_subtype)[,1]]
hc <- hclust(as.dist(samplePearsonDist), method="ward.D")
clusters=cutree(hc,k=4)

library(gplots)
heatmap = paste(output.folder, '/heatmap.png', sep="")
png(heatmap)
heatmap.2(samplePearsonDist, Rowv=as.dendrogram(hc),
          symm=TRUE, trace="none", col=colors,
          ColSideColors = cols,
          margins=c(1,5), labCol=FALSE )
dev.off()




#heatmap(euclidean)
sampleDists <- dist(t(logcounts(sce)[rownames(hvg.out)[1:100]]))
plot(hclust(sampleDists,method="ward.D")) 
library(RColorBrewer)
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(sce$molecular_subtype, sep="-" )
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
library("rafalib")
cols=palette(brewer.pal(8, "Dark2"))[as.data.frame(sce$molecular_subtype)[,1]]
hc <- hclust(sampleDists, method="ward.D")
library(gplots)
heatmap = paste(output.folder, '/heatmap.png', sep="")
png(heatmap)
heatmap.2( sampleDistMatrix, Rowv=as.dendrogram(hc),
           symm=TRUE, trace="none", col=colors,
           ColSideColors = cols,
           margins=c(1,5), labCol=FALSE )
dev.off()

#heatmap(spearman)
library(psych)
sampleSpearman=cor(logcounts(sce)[rownames(hvg.out)[1:100],],method="spearman")
sampleSpearmanDist=cor2dist(sampleSpearman)

rownames(sampleSpearmanDist) <- paste(sce$molecular_subtype, sep="-" )
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
plot(hclust(as.dist(sampleSpearmanDist), method="ward.D2"))
cols=palette(brewer.pal(8, "Dark2"))[as.data.frame(sce$molecular_subtype)[,1]]
hc <- hclust(as.dist(sampleSpearmanDist), method="ward.D")
library(gplots)
heatmap = paste(output.folder, '/heatmapSpearman.png', sep="")
png(heatmap)
heatmap.2(sampleSpearmanDist, Rowv=as.dendrogram(hc),
          symm=TRUE, trace="none", col=colors,
          ColSideColors = cols,
          margins=c(1,5), labCol=FALSE )
dev.off()

samplePearson=cor(logcounts(sce)[rownames(hvg.out)[1:100],],method="pearson")
samplePearsonDist=cor2dist(samplePearson)
plot(hclust(as.dist(samplePearsonDist),method="ward.D"))
rownames(samplePearsonDist) <- paste(sce$molecular_subtype, sep="-" )
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

cols=palette(brewer.pal(8, "Dark2"))[as.data.frame(sce$molecular_subtype)[,1]]
hc <- hclust(as.dist(samplePearsonDist), method="ward.D")
clusters=cutree(hc,k=4)

library(gplots)
heatmap = paste(output.folder, '/heatmap.png', sep="")
png(heatmap)
heatmap.2(samplePearsonDist, Rowv=as.dendrogram(hc),
          symm=TRUE, trace="none", col=colors,
          ColSideColors = cols,
          margins=c(1,5), labCol=FALSE )
dev.off()

#tsne
tsne = paste(output.folder, '/tsne.png', sep="")
png(tsne, width=1200, height=600)
out5 <- plotTSNE(sce, colour_by=top.hvg, run_args=list(perplexity=5,
                 feature_set=chosen)) + fontsize + ggtitle("Perplexity = 5")+ 
  coord_fixed(ratio=2 )
out10 <- plotTSNE(sce, colour_by=top.hvg, run_args=list(perplexity=10, 
                                     feature_set=chosen)) + fontsize + ggtitle("Perplexity = 10")+ 
  coord_fixed(ratio=2 )
out20 <- plotTSNE(sce, colour_by=top.hvg, run_args=list(perplexity=20, 
                                     feature_set=chosen)) + fontsize + ggtitle("Perplexity = 20")+ 
  coord_fixed(ratio=2 )
multiplot(out5, out10, out20, cols=3)
dev.off()

pca = paste(output.folder, '/pca.png', sep="")
png(pca)
pca1 <- plotPCA(sce, colour_by=top.hvg, run_args=list(exprs_values="norm_exprs")) + fontsize
multiplot(pca1)
dev.off()

#pca1 <- plotPCA(sce,exprs_values="exprs", colour_by=top.hvg) + fontsize
#multiplot(pca1)
#set.seed(100)
#plotPCA(sce, colour_by="total_features_by_counts") + fontsize
#saveRDS(file="hsc_data.rds", sce)
print("--PCA Complete")
chosen.exprs <- norm_exprs(sce)[chosen,]
my.dist <- dist(t(chosen.exprs))
my.tree <- hclust(my.dist, method="ward.D2")
library(dynamicTreeCut)
my.clusters <- unname(cutreeDynamic(my.tree, method="tree", verbose=0))


heat.vals <- chosen.exprs - rowMeans(chosen.exprs)
clust.col <- rainbow(max(my.clusters))
heatmap = paste(output.folder, '/heatmap.png', sep="")
png(heatmap)
heatmap.2(heat.vals, col=bluered, symbreak=TRUE, trace="none", cexRow=0.3, Colv=as.dendrogram(my.tree))
dev.off()
print("--Heatmap created")

library('Seurat')
markers2 <- FindMarkers(object = counts(dds), cells.1 = rownames(colData[colData[,specifiedFactor] == "2cell",,drop=F]),
                        cells.2 = rownames(colData[!(colData[,specifiedFactor] == "2cell"),,drop=F]), test.use = "DESeq2")
markers4 <- FindMarkers(object = counts(dds), cells.1 = rownames(colData[colData[,specifiedFactor] == "4cell",,drop=F]),
                        cells.2 = rownames(colData[!(colData[,specifiedFactor] == "4cell"),,drop=F]), test.use = "DESeq2")
markers8 <- FindMarkers(object = counts(dds), cells.1 = rownames(colData[colData[,specifiedFactor] == "8cell",,drop=F]),
                        cells.2 = rownames(colData[!(colData[,specifiedFactor] == "8cell"),,drop=F]), test.use = "DESeq2")
markers16 <- FindMarkers(object = counts(dds), cells.1 = rownames(colData[colData[,specifiedFactor] == "16cell",,drop=F]),
                         cells.2 = rownames(colData[!(colData[,specifiedFactor] == "16cell"),,drop=F]), test.use = "DESeq2")
markersblast <- FindMarkers(object = counts(dds), cells.1 = rownames(colData[colData[,specifiedFactor] == "blast",,drop=F]),
                            cells.2 = rownames(colData[!(colData[,specifiedFactor] == "blast"),,drop=F]), test.use = "DESeq2")