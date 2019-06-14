#' @importFrom methods setClass setMethod setValidity validObject new slot
NULL

#' The Startrac Class
#'
#' The Startrac object store the data for tcr-based T cell dynamics analyis. The slots contained in Startrac object are listed below:
#' @slot aid character. aid of the object, used for identification of the object. For example, patient id. default: "AID"
#' @slot cell.data data.frame. Each line for a cell, and these columns as required: `Cell_Name`, `clone.id`, `patient`, `majorCluster`, `loc`
#' @slot cell.perm.data object. list of `Startrac`` objects constructed from permutated cell data
#' @slot clonotype.data data.frame. Each line for a clonotype; contain the clonotype level indexes information
#' @slot cluster.data data.frame. Each line for a cluster; contain the cluster level indexes information
#' @slot pIndex.migr data.frame. Each line for a cluster; pairwise migration index between the two locations indicated in the column name.
#' @slot pIndex.tran data.frame. Each line for a cluster; pairwise transition index betwwen the two major clusters indicated by the row name and column name.
#' @slot cluster.sig.data data.frame. Each line for a cluster; contains the p values of cluster indices.
#' @slot pIndex.sig.migr data.frame. Each line for a cluster; contains the p values of pairwise migration indices.
#' @slot pIndex.sig.tran data.frame. Each line for a cluster; contains the p values of pairwise transition indices.
#' @slot clonotype.dist.loc matrix. Each line for a clonotype and describe the cells distribution among the locations.
#' @slot clonotype.dist.cluster matrix. Each line for a clonotype and describe the cells distribution among the clusters.
#' @slot clust.size array. Number of cells of each major cluster.
#' @slot patient.size array. Number of cells of each patient. 
#' @slot clone.size array. Number of cells of each clone.
#' @slot clone2patient array. Mapping from patient id to clone id.
#' @name Startrac
#' @rdname Startrac
#' @aliases Startrac-class
#' @exportClass Startrac
Startrac <- setClass("Startrac",
                     slots = c(aid = "character",
                               cell.data = "data.frame",
                               cell.perm.data = "list",
                               clonotype.data = "data.frame",
                               cluster.data = "data.frame",
                               cluster.sig.data = "data.frame",
                               pIndex.migr = "data.frame",
                               pIndex.tran = "data.frame",
                               pIndex.sig.migr = "data.frame",
                               pIndex.sig.tran = "data.frame",
                               clonotype.dist.loc = "matrix",
                               clonotype.dist.cluster = "matrix",
                               clust.size = "array",
                               patient.size = "array",
                               clone.size = "array",
                               clone2patient = "array"))


setValidity("Startrac",
            function(object) {
              msg <- NULL
              if(!is.data.frame(object@cell.data) || 
                 !all(c("Cell_Name","clone.id","patient","majorCluster","loc") %in% colnames(object@cell.data))){
                msg <- sprintf("cell.data must be data.frame and contain these columns: Cell_Name, clone.id, patient, majorCluster, loc")
              }
              if (is.null(msg)) TRUE
              else msg
            }
)


#' show method for Startrac
#
#' @param object A Startrac object
#' @name show
#' @aliases show,Startrac-method
#' @docType methods
#' @rdname show-methods
setMethod("show",
  signature = "Startrac",
  definition = function(object) {
    cat(sprintf("An object of class %s, aid %s, %d cells from %d clonotypes\n",
                class(object),
                object@aid,
                nrow(object@cell.data),
                nrow(object@clonotype.data)))
    invisible(x = NULL)
  }
)

#' initialize method for Startrac
#
#' @param .Object A Startrac object
#' @param cell.data data.frame contains the input data
#' @param aid character analysis id
#' @param n.perm integer number of permutation will be performed. If NULL, no permutation. (default: NULL)
#' @name initialize
#' @aliases initialize,Startrac-method
#' @docType methods
#' @rdname initialize-methods
#' @return an object of class \code{Startrac}
setMethod("initialize",
          signature = "Startrac",
          definition = function(.Object, cell.data, aid="AID",n.perm=NULL){
            .Object@aid <- aid
            .Object@cell.data <- as.data.frame(cell.data)
            validObject(.Object)
            .Object@clonotype.dist.loc <- unclass(table(.Object@cell.data[,c("clone.id","loc")]))
            .Object@clonotype.dist.cluster <- unclass(table(.Object@cell.data[,c("clone.id","majorCluster")]))
            
            .Object@clust.size <- unclass(table(.Object@cell.data$majorCluster))
            .Object@patient.size <- unclass(table(.Object@cell.data$patient))
            .Object@clone.size <- unclass(sort(table(.Object@cell.data$clone.id),decreasing = T))
            .clone2patient <- unique(.Object@cell.data[,c("patient","clone.id")])
            .Object@clone2patient <- as.array(structure(.clone2patient$patient, 
                                                        names=.clone2patient$clone.id))
            .Object@cell.perm.data <- list()
            if(!is.null(n.perm)){
              for(i in seq_len(n.perm)){
                perm.cell.data <- .Object@cell.data
                perm.cell.data$clone.id <- perm.cell.data$clone.id[sample(nrow(perm.cell.data))]
                .Object@cell.perm.data[[i]] <- new("Startrac",perm.cell.data,
                                                   aid=sprintf("perm%06d",i))
              }
            }
            return(.Object)
          }
)

#' Calculate cluster level indices
#
#' @name calIndex
#' @aliases calIndex calIndex,Startrac-method
#' 
#' @importFrom entropy entropy.empirical
#' @importFrom plyr llply
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @param object A Startrac object
#' @param n.perm integer number of permutation will be performed. If NULL, no permutation. (default: NULL)
#' @param cores number of core to be used. Passed to doParallel::registerDoParallel. default: NULL.
#' @return an object of class \code{Startrac}
Startrac.calIndex <- function(object,cores,n.perm) 
{
  ### cluster level expansion index (STARTRAC-expa)
  #### Todo: special case: number of clonotype is 1, i.e. sum(x>0)==1
  .entropy <- apply(object@clonotype.dist.cluster,2,function(x){ entropy.empirical(x,unit="log2") })
  .entropy.max <- apply(object@clonotype.dist.cluster,2,function(x){ log2(sum(x>0)) })
  object@cluster.data <- data.frame("aid"=object@aid,
                                    "majorCluster"=colnames(object@clonotype.dist.cluster),
                                    "expa"=1-.entropy/.entropy.max,
                                    stringsAsFactors = F)
  ### clone level migration and transition index
  object@clonotype.data <- data.frame("clone.id"=rownames(object@clonotype.dist.loc),
                                      "migr"=apply(object@clonotype.dist.loc,1,entropy.empirical,unit="log2"),
                                      "tran"=apply(object@clonotype.dist.cluster,1,entropy.empirical,unit="log2"))
  ### cluster level migration index (STARTRAC-migr) and transition index (STARTRAC-tran)
  object@cluster.data[["migr"]] <- apply(object@clonotype.dist.cluster,2,function(x){ sum(x*object@clonotype.data[names(x),"migr"])/sum(x) })
  object@cluster.data[["tran"]] <- apply(object@clonotype.dist.cluster,2,function(x){ sum(x*object@clonotype.data[names(x),"tran"])/sum(x) })
  if(!is.null(n.perm)){
    cl <- makeCluster(if(is.null(cores)) 1 else cores)
    registerDoParallel(cl)
    object@cell.perm.data <- llply(seq_len(n.perm),function(i){
      calIndex(object@cell.perm.data[[i]],cores=1)
    },.progress = "none",.parallel=T)
    stopCluster(cl)
    
  }
  return(object)
}

#' @export
setGeneric("calIndex", function(object,cores=NULL,n.perm=NULL) standardGeneric("calIndex"))

#' @rdname calIndex
#' @aliases calIndex
setMethod("calIndex", signature = "Startrac", definition = Startrac.calIndex)


#' Calculate pairwise indices
#
#' @name pIndex
#' @aliases pIndex pIndex,Startrac-method
#' 
#' @importFrom entropy entropy.empirical
#' @importFrom data.table dcast
#' @importFrom plyr ldply adply
#' @importFrom utils combn
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @param object A Startrac object
#' @param cores number of core to be used. Passed to doParallel::registerDoParallel. default: NULL.
#' @param n.perm integer number of permutation will be performed. If NULL, no permutation. (default: NULL)
#' @return an object of class \code{Startrac}
Startrac.pIndex <- function(object,cores,n.perm)
{
  cl <- makeCluster(if(is.null(cores)) 1 else cores)
  registerDoParallel(cl)
  ####### index given two cluster or loc
  ## migr 
  clone.dist.loc.majorCluster <- table(object@cell.data[,c("majorCluster","clone.id","loc")])
  ##tic("migr")
  withCallingHandlers({
    cls.migr.index.df <- ldply(seq_len(dim(clone.dist.loc.majorCluster)[1]),function(i,clone.dist.loc.majorCluster){
      dat.cls <- clone.dist.loc.majorCluster[i,,]
      i.name <- dimnames(clone.dist.loc.majorCluster)[["majorCluster"]][i]
      comb.loc <- as.data.frame(t(combn(colnames(dat.cls),2)),stringsAsFactors=F)
      dat.cls.pIndex.migr <- apply(comb.loc,1,function(x){
        dat.block <- dat.cls[,x]     
        dat.block.clone.index <- apply(dat.block,1,entropy::entropy.empirical,unit="log2")
        dat.block.clone.index[is.na(dat.block.clone.index)] <- 0
        sum(dat.block.clone.index*rowSums(dat.block)/sum(dat.block))
      })
      dat.cls.index <- cbind(data.frame(majorCluster=rep(i.name,nrow(comb.loc)),
                                     pIndex.migr=dat.cls.pIndex.migr,
                                     stringsAsFactors = F),
                          comb.loc)
      return(dat.cls.index)
    },clone.dist.loc.majorCluster=clone.dist.loc.majorCluster,.progress = "none",.parallel=T)
  },warning=function(w) {
    if(grepl("... may be used in an incorrect context:",conditionMessage(w)))
      ### strange bug, see https://github.com/hadley/plyr/issues/203
      invokeRestart("muffleWarning")
  })
  ##toc()
  cls.migr.index.df$crossLoc <- sprintf("%s-%s",cls.migr.index.df$V1,cls.migr.index.df$V2)
  object@pIndex.migr <- dcast(cls.migr.index.df,majorCluster ~ crossLoc,value.var = "pIndex.migr")
  object@pIndex.migr <- cbind(data.frame(aid=object@aid,stringsAsFactors = F),object@pIndex.migr)
  
  ## tran
  ##comb.cls <- as.data.frame(t(combn(colnames(object@clonotype.dist.cluster),2)),stringsAsFactors=F)
  comb.cls <- expand.grid(colnames(object@clonotype.dist.cluster),
                            colnames(object@clonotype.dist.cluster),stringsAsFactors = F)
  comb.cls <- comb.cls[comb.cls[,1]!=comb.cls[,2],]
  ##tic("tran")
  withCallingHandlers({
    cls.tran.index.df <- adply(comb.cls,1,function(x,object){
      dat.block <- object@clonotype.dist.cluster[,c(x[[1]],x[[2]])]
      dat.block.clone.index <- apply(dat.block,1,entropy::entropy.empirical,unit="log2")
      dat.block.clone.index[is.na(dat.block.clone.index)] <- 0
      data.frame(pIndex.tran=sum(dat.block.clone.index*rowSums(dat.block)/sum(dat.block)))
    },object=object,.progress = "none",.parallel=T)
  },warning=function(w) {
    if(grepl("... may be used in an incorrect context:",conditionMessage(w)))
      ### strange bug, see https://github.com/hadley/plyr/issues/203
      invokeRestart("muffleWarning")
  })
  stopCluster(cl)
  ##toc()
  
  object@pIndex.tran <- dcast(cls.tran.index.df,Var2~Var1,value.var = "pIndex.tran")
  colnames(object@pIndex.tran)[1] <- "majorCluster"
  object@pIndex.tran <- cbind(data.frame(aid=object@aid,stringsAsFactors = F),object@pIndex.tran)
  if(!is.null(n.perm)){
    for(i in seq_len(n.perm)){
      object@cell.perm.data[[i]] <- pIndex(object@cell.perm.data[[i]],cores=cores)
    }
  }
  return(object)  
}

#' @export
setGeneric("pIndex", function(object,cores=NULL,n.perm=NULL) standardGeneric("pIndex"))

#' @rdname pIndex
#' @aliases pIndex
setMethod("pIndex", signature = "Startrac", definition = Startrac.pIndex)




#' Get the p value given one Startrac object and a list of Startrac objects from permutation data
#
#' @name getSig
#' @aliases getSig getSig,Startrac-method
#' 
#' @importFrom plyr laply
#' @importFrom doParallel registerDoParallel
#' @importFrom data.table melt
#' @param obj A Startrac object
#' @param obj.perm A list of Startrac objects from permutation data 
#' @return an object of class \code{Startrac}
Startrac.getSig <- function(obj,obj.perm)
{
  .get.empirical.p <- function(obj,obj.perm,slot.name)
  {
    a.index <- melt(slot(obj,slot.name),id.vars=c("aid","majorCluster"),variable.name="index")
    perm.mtx <- t(laply(obj.perm,function(x){
      vv <- melt(slot(x,slot.name),id.vars=c("aid","majorCluster"),variable.name="index")
      vv.mtx <- as.matrix(vv[,"value",drop=F])
      rownames(vv.mtx) <- sprintf("%s.%s",vv[["majorCluster"]],vv[["index"]])
      colnames(vv.mtx) <- x@aid
      return(vv.mtx)
    }))
    stopifnot(all(sprintf("%s.%s",a.index$majorCluster,a.index$index)==rownames(perm.mtx)))
    a.index$p.value <- sapply(seq_along(a.index$value),function(i){
      v.rnd <- perm.mtx[i,]
      sum(v.rnd >= a.index$value[i])/length(v.rnd)
    })
    return(a.index)
  }
  obj@cluster.sig.data <- .get.empirical.p(obj,obj.perm,"cluster.data")
  obj@pIndex.sig.migr <- .get.empirical.p(obj,obj.perm,"pIndex.migr")
  obj@pIndex.sig.tran <- .get.empirical.p(obj,obj.perm,"pIndex.tran")
  return(obj)
}

#' @export
setGeneric("getSig", function(obj,obj.perm) standardGeneric("getSig"))

#' @rdname getSig
#' @aliases getSig
setMethod("getSig", signature = "Startrac", definition = Startrac.getSig)

