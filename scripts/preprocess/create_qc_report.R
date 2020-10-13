library(fastqcr)
args <- commandArgs(trailingOnly = TRUE)
qc.dir <-args[1]

qc <- qc_aggregate(qc.dir)
table <-summary(qc)

output = paste(qc.dir,"qc_summary.csv", sep="/")
write.csv(table , file =output)
