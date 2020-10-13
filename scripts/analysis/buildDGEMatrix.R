args <- commandArgs(trailingOnly = TRUE)

input.folder<- args[1] 
output.file<- args[2]
 
setwd(input.folder) 
filelist = list.files (pattern="*dge.txt.gz")
print(filelist)

#assuming tab separated values with a header
datalist = lapply(filelist, function(x)read.table(gzfile(x), header=TRUE, row.names=1, fill = FALSE))

#merge data files
dge_df=Reduce(function(x, y) merge(x, y, by='GENE', all=TRUE), datalist)

write.table(dge_df, file = gzfile(output.file), append = FALSE, quote = TRUE, sep="\t")