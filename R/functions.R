
#' dispaly message with time stamp
#' @param msg characters; message to display
#' @export
loginfo <- function(msg) {
  timestamp <- sprintf("%s", Sys.time())
  msg <- paste0("[",timestamp, "] ", msg,"\n")
  cat(msg)
}

#' warpper function for Startrac analysis
#' @importFrom entropy entropy.empirical
#' @importFrom data.table dcast
#' @importFrom plyr ldply adply llply
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom methods new slot
#' @param cell.data data.frame. Each line for a cell, and these columns as required: `Cell_Name`, `clone.id`, `patient`, `majorCluster`, `loc`
#' @param proj character. String used to annotate the project.
#' @param cores integer. number of core to be used. Passed to doParallel::registerDoParallel. default: NULL.
#' @param n.perm integer. number of permutation will be performed. If NULL, no permutation. (default: NULL)
#' @param verbose logical. wheter return intermediate result (some Startrac objects) 
#' @details run the Startrac pipeline
#' @return an list contains data.frame elements "cluster.data","pIndex.migr" and "pIndex.tran"
#' @export
#' @examples 
#' library("Startrac")
#' dat.file <- system.file("extdata/example.cloneDat.Zhang2018.txt",package = "Startrac")
#' in.dat <- read.table(dat.file,stringsAsFactors = FALSE,head=TRUE)
#' out <- Startrac.run(in.dat, proj="CRC", cores=2,verbose=FALSE)
#' 
Startrac.run <- function(cell.data, proj="CRC", cores=NULL,n.perm=NULL,verbose=F)
{
  ##tic("obj.proj")
  loginfo("initialize Startrac ...")
  obj.proj <- new("Startrac",cell.data,aid=proj,n.perm=n.perm)
  loginfo("calculate startrac index ...")
  obj.proj <- calIndex(obj.proj,cores=cores,n.perm=n.perm)
  loginfo("calculate pairwise index ...")
  obj.proj <- pIndex(obj.proj,cores=cores,n.perm=n.perm)
  if(!is.null(n.perm)){ 
    loginfo("get the significance")
    obj.proj <- getSig(obj.proj,obj.proj@cell.perm.data) 
  }
  ##toc()
  
  obj.list <- NULL
  if(length(obj.proj@patient.size)>1)
  {
    loginfo("calculate indices of each patient ...")
    patient.vec <- names(obj.proj@patient.size[obj.proj@patient.size > 30])
    cl <- makeCluster(if(is.null(cores)) 1 else cores)
    registerDoParallel(cl)
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
    stopCluster(cl)
    ##toc()
  }
  loginfo("collect result")
  ret <- list()
  ## cluster index
  ret.slot.names <- c("cluster.data","pIndex.migr","pIndex.tran")
  if(!is.null(n.perm)){ 
    ret.slot.names <- c(ret.slot.names,
                        c("cluster.sig.data","pIndex.sig.migr","pIndex.sig.tran")) 
  }
  for(v in ret.slot.names)
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
  loginfo("return")
  return(ret)
}
