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
filelist = list.files(pattern="*dge.txt.gz")

metadf <- read.table(metadata.file, header = TRUE, sep = "\t")

#rename column name with source ID
rename_columns<-function(file) {
  dataframe=read.table(gzfile(file), header=T, sep="\t")
  column_list=c('GENE') 
  runnum = gsub("_dge.txt.gz","",file)
  sourcedf <- fn$sqldf("select ID from metadf where Run='$runnum' ")
  sourcename <- as.character((sourcedf[1,1])) 
  for (i in (2:ncol(dataframe))) {
    column_name = paste(sourcename, i-1, sep="_")
    column_list=c(column_list,column_name)
  }

  #assign new column name to dataframe
  colnames(dataframe)=column_list   
  return (dataframe)
}

datalist = lapply(filelist, function(x){
  rename_columns(x)
})

#merge data files
dge_df=Reduce(function(x, y) merge(x, y, by='GENE', all=TRUE), datalist)

# replace na with 0
dge_df[is.na(dge_df)]=0

output.file = paste(output.folder, output.file.name ,sep="/" )
write.table(dge_df, gzfile(output.file), append = FALSE, quote = TRUE, sep="\t")

#write.table(dgeframe, file = output.file, append = FALSE, quote = TRUE, sep="\t")
