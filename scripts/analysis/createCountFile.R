#create count matrix for each individual

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

#assuming tab separated values without a header
datalist = lapply(filelist, function(x){
 tmp <- try(read.table( x, header=F, row.names=1 ))
  if (!inherits(tmp, 'try-error')) tmp
})

#assuming the same header/columns for all files
datafr = do.call("cbind", datalist)

#define an empty column_list vector
column_list=c()

#read metadata to find the ID based on the run number
metadf <- read.table(metadata.file, header = TRUE, sep = "\t")
for (i in (1:ncol(datafr))){
  run_numbers=as.vector(cbind(gsub("_results.txt","",gsub("htseq_","",filelist))))
  id=run_numbers[i]
  print (id)
  sourcedf <- fn$sqldf("select ID from metadf where Run='$id' ")
  sourcename <- as.character((sourcedf[1,1]))
  column_name=sourcename

  #combine column name into column_list vector
  column_list=c(column_list,column_name)
}

#assign column name to datafr
colnames(datafr)=column_list
print(column_list)
dim(datafr)
 
output.file = paste(output.folder, output.file.name ,sep="/" )
write.table(datafr, gzfile(output.file), append = FALSE, quote = TRUE, sep="\t")
             
