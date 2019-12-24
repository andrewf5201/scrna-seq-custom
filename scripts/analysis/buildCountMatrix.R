args <- commandArgs(trailingOnly = TRUE)

input.folder<- args[1] 
output.file<- args[2]


# Set the working directory
setwd(input.folder)
filelist = list.files(pattern="*counts.txt")
print(filelist)

#assuming tab separated values with a header
datalist = lapply(filelist, function(x)read.table( x, header=TRUE, row.names=1, fill = FALSE))
cts = do.call("cbind", datalist)
#remove last 5 rows which contains summary info
cts <- cts[1:(nrow(cts)-5),]

write.table(cts, file = output.file, append = FALSE, quote = TRUE, sep="\t")

