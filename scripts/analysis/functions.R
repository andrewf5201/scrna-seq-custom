library('properties')
get_dir_info <- function (algorithm, config, packagename) { 
  pipeline.dir = config$pipeline_dir
  experiment = config$experiment 
  
  analysis.dir = paste(pipeline.dir, "analysis", experiment, sep =
                         "/")
  input.dir = paste(analysis.dir, "data", sep = "/")
  output.dir.main = paste(analysis.dir, "results", packagename, sep =
                            "/")
  subdir = paste(Sys.Date(), algorithm, sep = '/')
  output.dir = paste(output.dir.main, subdir, sep = "/")
  dir.create(
    file.path(output.dir.main, subdir),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  dirList <- list("input.dir" = input.dir, "out.dir" = output.dir)
  return (dirList)
}

get_factors_info <-function (metadata.file) {
  df <- read.table(metadata.file, header = TRUE, sep = "\t")
  #get all factors defined in metadata file (from the 6th column)
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
  factor_info <-list("factorlists"=factorlists, "factors"=fctrs)
  return (factor_info)
} 

convertIDs <- function( ids, from, to, db, ifMultiple=c("putNA", "useFirst")) {
  stopifnot( inherits( db, "AnnotationDb" ) )
  ifMultiple <- match.arg( ifMultiple )
  suppressWarnings( selRes <- AnnotationDbi::select(
   db, keys=ids, keytype=from, columns=c(from,to) ) )
  if ( ifMultiple == "putNA" ) {
   duplicatedIds <- selRes[ duplicated( selRes[,1] ), 1 ]
   selRes <- selRes[ ! selRes[,1] %in% duplicatedIds, ]
  }
  return( selRes[ match( ids, selRes[,1] ), 2 ] )
}