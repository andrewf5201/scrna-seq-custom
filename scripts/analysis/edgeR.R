library(edgeR) 

args <- commandArgs(trailingOnly = TRUE)

algorithm<- args[1]
input.folder<- args[2]
output.folder<- args[3]
#convert to lowercase
algorithm <- tolower(algorithm)

# read input files
setwd(input.folder) 
filelist = list.files(pattern="*counts.txt")

#assuming tab separated values with a header   
datalist = lapply(filelist, function(x)read.table( x, header=TRUE, row.names=1, fill = FALSE)) 
cts = do.call("cbind", datalist) 
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

group= interaction(factorlists)
dgList <- DGEList(counts=cts, genes=rownames(cts), group=group) 

#Quality Assessment
pseudoCounts <-log2(dgList$counts +1)
head(pseudoCounts)
hist(pseudoCounts) 
 
dgList[is.na(dgList)] <- 0 
apply(dgList$counts,2,function(x) is.na(x))  

#filtering
countsPerMillion <- cpm(dgList)
summary(countsPerMillion)

countCheck <- countsPerMillion > 1
head(countCheck)
keep <- which(rowSums(countCheck) >= 2)
dgList <- dgList[keep,]
summary(cpm(dgList))

dgList$samples$lib.size
dgList$samples$lib.size <- colSums(dgList$counts)
dgList$samples
dgList$samples$lib.size

#normalization
if (algorithm == "default") { 
  #default normalization
  dgList <- calcNormFactors(dgList)
  nc <- cpm(dgList, normalized.lib.sizes=FALSE)
  
  defaultNormPlotName=paste(output.folder, '/edgeR/norm_default.png', sep="")
  png(defaultNormPlotName)
  plotMDS.DGEList(nc, col = as.numeric(group))
  dev.off()
} else if (algorithm == "tmm") { 
  #tmm normalization
  dgList <- calcNormFactors(dgList, method="TMM")
  
  tmmPlotName = paste(output.folder, '/edgeR/norm_tmm.png', sep="")
  png(tmmPlotName)
  plotMDS.DGEList(dgList)
  dev.off()
} else if (algorithm == "cpm") { 
  #cpm normalization 
  cpm=cpm(dgList)
  cpm=cpm(dgList,lib.size=dgList$samples$lib.size,log = TRUE)
  
  cpmPlotName = paste(output.folder, '/edgeR/norm_cpm.png', sep="")
  png(cpmPlotName) 
  plotMDS.DGEList(cpm)
  dev.off()
} else if (algorithm == "rpkm") {  
  gene_len <-read.delim(paste(input.folder,"hg38_gene_length.txt",sep=""), header=F, sep="\t")
  colnames(gene_len) <- c("GeneName","Len") 
  m <- match(rownames(dgList), gene_len$GeneName)
  gene.lengths <- gene_len$Len[m]
  
  #rpkm normalization
  dgList<- rpkm(dgList, gene.lengths)  
  
  rpkmPlotName = paste(output.folder, 'edgeR/norm_rpkm.png', sep="")
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
meanPlotName = paste(output.folder, '/edgeR/meanVar.png', sep="")
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
bcvPlot = paste(output.folder, '/edgeR/bcv.png', sep="")
png(bcvPlot) 
plotBCV(dgList)  
dev.off()
 
fit <- glmFit(dgList)
lrt <- glmLRT(fit, coef=4)
edgeR_result <- topTags(lrt)
deGenes <- decideTestsDGE(lrt, p=0.001)
deGenes <- rownames(lrt)[as.logical(deGenes)]
plotSmear(lrt, de.tags=deGenes)
abline(h=c(-1, 1), col=2)



