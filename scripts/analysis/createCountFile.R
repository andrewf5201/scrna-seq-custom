install.packages("sqldf",repos = "http://cran.us.r-project.org")

library(sqldf)
args <- commandArgs(trailingOnly = TRUE)

input.folder<- args[1] 
output.folder<- args[2]
output.file.name<- args[3] 
metadata.file<- args[4]

 
# Set the working directory
setwd(input.folder)
print(input.folder)
filelist = list.files(pattern="*.txt")

#assuming tab separated values with a header    
datalist = lapply(filelist, function(x){
 tmp <- try(read.table( x, header=F, row.names=1 ))
  if (!inherits(tmp, 'try-error')) tmp
})

#assuming the same header/columns for all files
datafr = do.call("cbind", datalist)
column_list=c()
df <- read.table(metadata.file, header = TRUE, sep = "\t")
for (i in (1:ncol(datafr))){
  id_number=as.vector(cbind(gsub("_results.txt","",gsub("htseq_","",filelist))))
  id=id_number[i] 
  sourcedf <- fn$sqldf("select ID from df where Run='$id' ")
  sourcename <- as.character((sourcedf[1,1]))
  column_name=sourcename
  column_list=c(column_list,column_name)
}
colnames(datafr)=column_list
print(column_list)
dim(datafr)
 
output.file = paste(output.folder, output.file.name ,sep="/" )
write.table(datafr, file = output.file, append = FALSE, quote = TRUE, sep="\t")
             
