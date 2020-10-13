library(edgeR) 

args <- commandArgs(trailingOnly = TRUE)

#get input values
args <- commandArgs(trailingOnly = TRUE)
algorithm <-args[1]
config_file <- args[2]
script_dir <-args[3]

#convert to lowercase
algorithm <- tolower(algorithm)

function_script = paste(script_dir, "functions.R", sep='/')
source(function_script)

config <- read.properties(config_file)
dirlist = get_dir_info(algorithm, config, packagename="edgeR")
input.folder = dirlist$input.dir
output.folder= dirlist$out.dir
count.matrix.file = config$matrix_file
metadata.file = config$metadata_file
gene.length.file = config$gene_length_file
specifiedFactor = config$cluster_factor

sprintf("--loading count matirx file %s",count.matrix.file)

if (endsWith(count.matrix.file, "gz")) {
  cts <- read.table(gzfile(count.matrix.file), sep="\t", header=TRUE)
} else {
  cts <- read.table(count.matrix.file, sep="\t", header=TRUE)
}
print(cts)

factor_info=get_factors_info(metadata.file)
fctrs=factor_info$factors
factorlists=factor_info$factorlists
colData=as.data.frame(factorlists)
colnames(colData)=fctrs
rownames(colData)=colnames(cts)
print("colData created")

#--create DGEList
group= interaction(factorlists)
dgList <- DGEList(counts=cts, genes=rownames(cts), group=group) 

#Quality Assessment
pseudoCounts <-log2(dgList$counts +1)
head(pseudoCounts)

hist(pseudoCounts) 

dgList[is.na(dgList)] <- 0 
apply(dgList$counts,2,function(x) is.na(x))  

#gene filtering by CPM
countsPerMillion <- cpm(dgList)


summary.file = paste(output.folder, "cpmSum1.txt" ,sep="/" )
sink(summary.file)
summary(countsPerMillion)
sink()


countCheck <- countsPerMillion > 1
head(countCheck)
keep <- which(rowSums(countCheck) >= 2)
dgList <- dgList[keep,]

#write summary to a file
summary.file = paste(output.folder, "cpmSumdgList.txt" ,sep="/" )
sink(summary.file)
summary(cpm(dgList))
sink()


dgList$samples$lib.size
dgList$samples$lib.size <- colSums(dgList$counts)
dgList$samples
dgList$samples$lib.size

#normalization
if (algorithm == "default") { 
  #default normalization
  dgList <- calcNormFactors(dgList)
  nc <- cpm(dgList, normalized.lib.sizes=FALSE)
  
  defaultNormPlotName=paste(output.folder, '/norm_default.png', sep="")
  png(defaultNormPlotName)
  plotMDS.DGEList(nc, col = as.numeric(group))
  dev.off()
} else if (algorithm == "tmm") { 
  #tmm normalization
  dgList <- calcNormFactors(dgList, method="TMM")
  
  tmmPlotName = paste(output.folder, '/norm_tmm.png', sep="")
  png(tmmPlotName)
  plotMDS.DGEList(dgList)
  dev.off()
} else if (algorithm == "cpm") { 
  #cpm normalization 
  cpm=cpm(dgList)
  cpm=cpm(dgList,lib.size=dgList$samples$lib.size,log = TRUE)
  
  cpmPlotName = paste(output.folder, '/norm_cpm.png', sep="")
  png(cpmPlotName) 
  plotMDS.DGEList(cpm)
  dev.off()
} else if (algorithm == "rpkm") {
  gene_len <- read.delim(gene.length.file, header=F, sep="\t")
  colnames(gene_len) <- c("GeneName","Len") 
  m <- match(rownames(dgList), gene_len$GeneName)
  gene.lengths <- gene_len$Len[m]
  
  #rpkm normalization
  dgList<- rpkm(dgList, gene.lengths)  
  
  rpkmPlotName = paste(output.folder, '/norm_rpkm.png', sep="")
  png(rpkmPlotName)
  plotMDS.DGEList(dgList)
  dev.off()
}

dgList <- estimateCommonDisp(dgList)
names(dgList)
dgList$common.dispersion

dgList <- estimateTagwiseDisp( dgList , prior.n = 10 )
summary( dgList$tagwise.dispersion )

#meanVarPlot
meanPlotName = paste(output.folder, '/meanVar.png', sep="")
png(meanPlotName)
meanVarPlot <- plotMeanVar( dgList , show.raw.vars=TRUE ,
                            show.tagwise.vars=TRUE ,
                            show.binned.common.disp.vars=FALSE ,
                            show.ave.raw.vars=FALSE ,
                            dispersion.method = "qcml" , NBline = TRUE ,
                            nbins = 100 ,
                            pch = 16 ,
                            xlab ="Mean Expression (Log10 Scale)" ,
                            ylab = "Variance (Log10 Scale)" ,
                            main = "Mean-Variance Plot" )
dev.off()

#bcvPlot
bcvPlot = paste(output.folder, '/bcv.png', sep="")
png(bcvPlot) 
plotBCV(dgList)  
dev.off()

fit <- glmFit(dgList)
lrt <- glmLRT(fit, coef=1)

#TODO Save result to file??
edgeR_result <- topTags(lrt)
deGenes <- decideTestsDGE(lrt, p=0.001)
deGenes <- rownames(lrt)[as.logical(deGenes)]
#write summary to a file
DElist = paste(output.folder, "deGenes.txt" ,sep="/" )
sink(DElist)
deGenes
sink()

#TODO: save to file
smear = paste(output.folder, '/plotSmear.png', sep="")
png(smear) 
plotSmear(lrt, de.tags=deGenes)
abline(h=c(-1, 1), col=2)
dev.off



