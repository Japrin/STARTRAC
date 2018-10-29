#' warpper function for Startrac analysis
#' @importFrom entropy entropy.empirical
#' @importFrom data.table dcast
#' @importFrom plyr ldply adply llply
#' @importFrom doParallel registerDoParallel
#' @importFrom methods new slot
#' @param cell.data data.frame. Each line for a cell, and these columns as required: `Cell_Name`, `clone.id`, `patient`, `majorCluster`, `loc`
#' @param proj character. String used to annotate the project.
#' @param cores integer. number of core to be used. Passed to doParallel::registerDoParallel. default: NULL.
#' @param verbose logical. wheter return intermediate result (some Startrac objects) 
#' @details run the Startrac pipeline
#' @return an list contains data.frame elements "cluster.data","pIndex.migr" and "pIndex.tran"
#' @export
#' @examples 
#' library("Startrac")
#' dat.file <- system.file("extdata/example.cloneDat.Zhang2018.txt",package = "Startrac")
#' in.dat <- read.table(dat.file,stringsAsFactors = F,head=T)
#' out <- Startrac.run(in.dat, proj="CRC", cores=NULL,verbose=F)
#' 
Startrac.run <- function(cell.data, proj="CRC", cores=NULL,verbose=F)
{
  ##tic("obj.proj")
  obj.proj <- new("Startrac",cell.data,aid=proj)
  obj.proj <- calIndex(obj.proj)
  obj.proj <- pIndex(obj.proj,cores=cores)
  ##toc()
  
  obj.list <- NULL
  if(length(obj.proj@patient.size)>1)
  {
    patient.vec <- names(obj.proj@patient.size[obj.proj@patient.size > 30])
    registerDoParallel(cores = cores)
    ##tic("obj.list")
    withCallingHandlers({
      obj.list <- llply(patient.vec,function(pid,cell.data){
        require("Startrac")
        obj <- new("Startrac",subset(cell.data,patient==pid),aid=pid)
        obj <- calIndex(obj)
        obj <- pIndex(obj,cores=1)
        obj
      },cell.data=cell.data,.progress = "none",.parallel=T)
    },warning=function(w) {
      if(grepl("... may be used in an incorrect context:",conditionMessage(w)))
        ### strange bug, see https://github.com/hadley/plyr/issues/203
        invokeRestart("muffleWarning")
    })
    ##toc()
  }
  ret <- list()
  ## cluster index
  for(v in c("cluster.data","pIndex.migr","pIndex.tran"))
  {
    ret[[v]] <- slot(obj.proj,v)
    if(!is.null(obj.list)){
      ret[[v]] <- rbind(ret[[v]],ldply(obj.list,function(obj){
        slot(obj,v)
      }))
    }
  }
  if(verbose){
    ret[["objects"]] <- c(obj.proj,obj.list)
  }
  return(ret)
}
