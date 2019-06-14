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
#' @importFrom plyr llply
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom doParallel registerDoParallel
#' @param object A Startrac object
#' @param n.perm integer number of permutation will be performed. If NULL, no permutation. (default: NULL)
#' @param cores number of core to be used. Passed to doParallel::registerDoParallel. default: NULL.
#' @param normEntropy logical; whether normalize migration and transition index. default: FALSE.
#' @return an object of class \code{Startrac}
Startrac.calIndex <- function(object,cores,n.perm,normEntropy)
{
  ### cluster level expansion index (STARTRAC-expa)
  #### Todo: special case: number of clonotype is 1, i.e. sum(x>0)==1
  .entropy <- mcol.entropy(object@clonotype.dist.cluster)
  .entropy.max <- log2(colSums(object@clonotype.dist.cluster > 0))
  object@cluster.data <- data.frame("aid"=object@aid,
                                    "majorCluster"=colnames(object@clonotype.dist.cluster),
                                    "expa"=1-.entropy/.entropy.max,
                                    stringsAsFactors = F)
  ### clone level migration and transition index
  if(normEntropy){
    .entropy.migr.max <- log2(ncol(object@clonotype.dist.loc))
    .entropy.tran.max <- log2(ncol(object@clonotype.dist.cluster))
  }else{
    .entropy.migr.max <- 1
    .entropy.tran.max <- 1
  }
  object@clonotype.data <- data.frame("clone.id"=rownames(object@clonotype.dist.loc),
                                      "migr"=mrow.entropy(object@clonotype.dist.loc)/.entropy.migr.max,
                                      "tran"=mrow.entropy(object@clonotype.dist.cluster)/.entropy.tran.max)
  ### cluster level migration index (STARTRAC-migr) and transition index (STARTRAC-tran)
  weights.mtx <- sweep(object@clonotype.dist.cluster,2,colSums(object@clonotype.dist.cluster),"/")
  index.mtx <- t(weights.mtx) %*% (as.matrix(object@clonotype.data[,c("migr","tran")]))
  object@cluster.data <- cbind(object@cluster.data,index.mtx)
  if(!is.null(n.perm)){
    #cl <- makeCluster(if(is.null(cores)) (detectCores()-2) else cores)
    #registerDoParallel(cl)
    registerDoParallel(if(is.null(cores)) (detectCores()-2) else cores)
    object@cell.perm.data <- llply(object@cell.perm.data,function(x){
      calIndex(x,cores=1,normEntropy=normEntropy)
    },.progress = "none",.parallel=T)
    #stopCluster(cl)
    
  }
  return(object)
}

#' @export
setGeneric("calIndex", function(object,cores=NULL,n.perm=NULL,normEntropy=FALSE) standardGeneric("calIndex"))

#' @rdname calIndex
#' @aliases calIndex
setMethod("calIndex", signature = "Startrac", definition = Startrac.calIndex)


#' Calculate pairwise indices
#
#' @name pIndex
#' @aliases pIndex pIndex,Startrac-method
#'
#' @importFrom data.table dcast
#' @importFrom plyr ldply adply
#' @importFrom utils combn
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom doParallel registerDoParallel
#' @param object A Startrac object
#' @param cores number of core to be used. Passed to doParallel::registerDoParallel. default: NULL.
#' @param n.perm integer number of permutation will be performed. If NULL, no permutation. (default: NULL)
#' @return an object of class \code{Startrac}
Startrac.pIndex <- function(object,cores,n.perm)
{
  ####### index given two cluster or loc
  ## migr 
  clone.dist.loc.majorCluster <- table(object@cell.data[,c("majorCluster","clone.id","loc")])
  
  #tic("migr")
  withCallingHandlers({
    cls.migr.index.df <- ldply(seq_len(dim(clone.dist.loc.majorCluster)[1]),function(i,clone.dist.loc.majorCluster){
      dat.cls <- clone.dist.loc.majorCluster[i,,]
      i.name <- dimnames(clone.dist.loc.majorCluster)[["majorCluster"]][i]
      comb.loc <- as.data.frame(t(combn(colnames(dat.cls),2)),stringsAsFactors=F)
      dat.cls.pIndex.migr <- apply(comb.loc,1,function(x){
        dat.block <- dat.cls[,x]
        dat.block.clone.index <- mrow.entropy(dat.block)
        dat.block.clone.index[is.na(dat.block.clone.index)] <- 0
        t(rowSums(dat.block)/sum(dat.block)) %*% dat.block.clone.index
      })
      dat.cls.index <- cbind(data.frame(majorCluster=rep(i.name,nrow(comb.loc)),
                                     pIndex.migr=dat.cls.pIndex.migr,
                                     stringsAsFactors = F),
                          comb.loc)
      return(dat.cls.index)
    },clone.dist.loc.majorCluster=clone.dist.loc.majorCluster,.progress = "none",.parallel=F)
  },warning=function(w) {
    if(grepl("... may be used in an incorrect context:",conditionMessage(w)))
      ### strange bug, see https://github.com/hadley/plyr/issues/203
      invokeRestart("muffleWarning")
  })
  #toc()
  
  cls.migr.index.df$crossLoc <- sprintf("%s-%s",cls.migr.index.df$V1,cls.migr.index.df$V2)
  object@pIndex.migr <- dcast(cls.migr.index.df,majorCluster ~ crossLoc,value.var = "pIndex.migr")
  object@pIndex.migr <- cbind(data.frame(aid=object@aid,stringsAsFactors = F),object@pIndex.migr)
  
  ## tran
  comb.cls <- expand.grid(colnames(object@clonotype.dist.cluster),
                            colnames(object@clonotype.dist.cluster),stringsAsFactors = F)
  comb.cls <- comb.cls[comb.cls[,1]!=comb.cls[,2],]
  
  #tic("tran")
  withCallingHandlers({
    cls.tran.index.df <- adply(comb.cls,1,function(x,object){
      dat.block <- object@clonotype.dist.cluster[,c(x[[1]],x[[2]])]
      dat.block.clone.index <- mrow.entropy(dat.block)
      dat.block.clone.index[is.na(dat.block.clone.index)] <- 0
      data.frame(pIndex.tran= t(rowSums(dat.block)/sum(dat.block)) %*% dat.block.clone.index)
    },object=object,.progress = "none",.parallel=F)
  },warning=function(w) {
    if(grepl("... may be used in an incorrect context:",conditionMessage(w)))
      ### strange bug, see https://github.com/hadley/plyr/issues/203
      invokeRestart("muffleWarning")
  })
  #toc()
  
  object@pIndex.tran <- dcast(cls.tran.index.df,Var2~Var1,value.var = "pIndex.tran")
  colnames(object@pIndex.tran)[1] <- "majorCluster"
  object@pIndex.tran <- cbind(data.frame(aid=object@aid,stringsAsFactors = F),object@pIndex.tran)
  if(!is.null(n.perm)){
    #cl <- makeCluster(if(is.null(cores)) (detectCores()-2)  else cores)
    #registerDoParallel(cl)
    registerDoParallel(if(is.null(cores)) (detectCores()-2)  else cores)
    object@cell.perm.data <- llply(object@cell.perm.data,function(x){
      pIndex(x,n.perm=NULL)
    },.progress = "none",.parallel=T)
    #stopCluster(cl)
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
    if(!is.null(obj.perm)){
      perm.mtx <- t(laply(obj.perm,function(x){
        vv <- melt(slot(x,slot.name),id.vars=c("aid","majorCluster"),variable.name="index")
        vv.mtx <- as.matrix(vv[,"value",drop=F])
        rownames(vv.mtx) <- sprintf("%s.%s",vv[["majorCluster"]],vv[["index"]])
        colnames(vv.mtx) <- x@aid
        return(vv.mtx)
      }))
      stopifnot(all(sprintf("%s.%s",a.index$majorCluster,a.index$index)==rownames(perm.mtx)))
      a.index$p.value <- sapply(seq_along(a.index$value),function(i){
        v.rnd <- perm.mtx[i,,drop=F]
        sum(v.rnd >= a.index$value[i])/length(v.rnd)
      })
    }else{
      a.index$p.value <- NA
    }
    return(a.index)
  }
  obj@cluster.sig.data <- .get.empirical.p(obj,obj.perm,"cluster.data")
  obj@pIndex.sig.migr <- .get.empirical.p(obj,obj.perm,"pIndex.migr")
  obj@pIndex.sig.tran <- .get.empirical.p(obj,obj.perm,"pIndex.tran")
  return(obj)
}

#' @export
setGeneric("getSig", function(obj,obj.perm=NULL) standardGeneric("getSig"))

#' @rdname getSig
#' @aliases getSig
setMethod("getSig", signature = "Startrac", definition = Startrac.getSig)



#' The StartracOUt Class
#'
#' Object store the result of Startrac.run:
#' @slot proj character. identification of the object. For example, patient id. default: "AID"
#' @slot cluster.data data.frame. Each line for a cluster; contain the cluster level indexes information
#' @slot pIndex.migr data.frame. Each line for a cluster; pairwise migration index between the two locations indicated in the column name.
#' @slot pIndex.tran data.frame. Each line for a cluster; pairwise transition index betwwen the two major clusters indicated by the row name and column name.
#' @slot cluster.sig.data data.frame. Each line for a cluster; contains the p values of cluster indices.
#' @slot pIndex.sig.migr data.frame. Each line for a cluster; contains the p values of pairwise migration indices.
#' @slot pIndex.sig.tran data.frame. Each line for a cluster; contains the p values of pairwise transition indices.
#' @slot objects list. other objects
#' @name StartracOut
#' @rdname StartracOut
#' @aliases StartracOut-class
#' @exportClass StartracOut
StartracOut <- setClass("StartracOut",
                     slots = c(proj = "character",
                               cluster.data = "data.frame",
                               cluster.sig.data = "data.frame",
                               pIndex.migr = "data.frame",
                               pIndex.tran = "data.frame",
                               pIndex.sig.migr = "data.frame",
                               pIndex.sig.tran = "data.frame",
                               objects = "list"))

#' initialize method for StartracOut
#
#' @param .Object A StartracOut object
#' @param proj character analysis id
#' @aliases initialize,StartracOut-method
#' @docType methods
#' @return an object of class \code{StartracOut}
setMethod("initialize",
          signature = "StartracOut",
          definition = function(.Object, proj="AID"){
            .Object@proj <- proj
            return(.Object)
          }
)


#' show method for StartracOut
#
#' @importFrom utils head
#' @param object A StartracOut object
#' @aliases show,StartracOut-method
#' @docType methods
setMethod("show",
          signature = "StartracOut",
          definition = function(object) {
            cat(sprintf("An object of class %s, proj %s:\n",
                        class(object),
                        object@proj))
            cat("head of the clusters' index:\n")
            print(head(object@cluster.data))
            cat("head of the pairwise migration index:\n")
            print(head(object@pIndex.migr))
            cat("head of the pairwise transition index:\n")
            print(head(object@pIndex.tran[,1:min(5,ncol(object@pIndex.tran))]))
            invisible(x = NULL)
          }
)

#' plot the indexes
#
#' @name plot
#' @aliases plot plot,StartracOut-method
#' 
#' @importFrom plyr laply
#' @importFrom doParallel registerDoParallel
#' @importFrom data.table melt as.data.table melt
#' @importFrom ggpubr ggbarplot ggboxplot
#' @importFrom ggplot2 facet_wrap theme element_text aes geom_text
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom RColorBrewer brewer.pal
#' @importFrom circlize colorRamp2
#' @importFrom grDevices colorRampPalette
#' @param obj A object of StartracOut
#' @param index.type one of "cluster.all", "pairwise.migr", "pairwise.tran". (default:"cluster.all")
#' @param byPatient logical. plot indexes of each patient (default: FALSE)
#' @return a ggplot2 object or Heatmap-class object
StartracOut.plot <- function(obj,index.type,byPatient)
{
  if(index.type=="cluster.all"){
    if(byPatient){
      p <- ggboxplot(as.data.table(obj@cluster.sig.data)[aid!=obj@proj,],
                x="majorCluster",y="value",palette = "npg",
                color = "index", add = "point", outlier.colour=NULL) +
        facet_wrap(~index,ncol=1,scales = "free_y") +
        theme(axis.text.x=element_text(angle = 60,hjust = 1))
    }else{
      dat.plot <- as.data.table(obj@cluster.sig.data)[aid==obj@proj,]
      dat.plot$p.value.label <- ""
      dat.plot$p.value.label[dat.plot$p.value < 0.05] <- "*"
      dat.plot$p.value.label[dat.plot$p.value < 0.01] <- "**"
      dat.plot$p.value.label[dat.plot$p.value < 0.001] <- "***"
      p <- ggbarplot(dat.plot,
                    x="majorCluster",y="value",palette = "npg",fill = "index") +
        facet_wrap(~index,ncol=1,scales = "free_y") +
        theme(axis.text.x=element_text(angle = 60,hjust = 1))
      if(!all(is.na(dat.plot$p.value))){
        p <- p + geom_text(aes(label=p.value.label,y=value+0.01),size=5)
      }
    }

  }else if(index.type=="pairwise.migr"){
    if(byPatient){
      p <- ggboxplot(as.data.table(obj@pIndex.sig.migr)[aid!=obj@proj,],
                     x="majorCluster",y="value",palette = "npg",
                     color = "index", add = "point", outlier.colour=NULL) +
        facet_wrap(~index,ncol=1,scales = "free_y") +
        theme(axis.text.x=element_text(angle = 60,hjust = 1))      
    }else{
      dat.plot <- as.data.table(obj@pIndex.sig.migr)[aid==obj@proj,]
      dat.plot$p.value.label <- ""
      dat.plot$p.value.label[dat.plot$p.value < 0.05] <- "*"
      dat.plot$p.value.label[dat.plot$p.value < 0.01] <- "**"
      dat.plot$p.value.label[dat.plot$p.value < 0.001] <- "***"
      p <- ggbarplot(dat.plot,
                x="majorCluster",y="value",palette = "npg",fill = "index") +
        facet_wrap(~index,ncol=1,scales = "free_y") +
        theme(axis.text.x=element_text(angle = 60,hjust = 1))
      if(!all(is.na(dat.plot$p.value))){
        p <- p + geom_text(aes(label=p.value.label,y=value+0.01),size=5)
      }
    }
  }else if(index.type=="pairwise.tran"){
    dat.plot <- as.matrix(subset(obj@pIndex.tran,aid==obj@proj)[,c(-1,-2)])
    rownames(dat.plot) <- subset(obj@pIndex.tran,aid==obj@proj)[,2]
    dat.plot[is.na(dat.plot)] <- 0
    col.heat <- colorRamp2(seq(0,0.12,length=15),
                           colorRampPalette(rev(brewer.pal(n=7,name="RdBu")))(15),
                           space = "LAB")
    p <- Heatmap(dat.plot,name="pIndex.tran",col = col.heat)
  }
  return(p)
}

#' @export
setGeneric("plot", function(obj,index.type="cluster.all",byPatient=F) standardGeneric("plot"))

#' @rdname plot
#' @aliases plot
setMethod("plot", signature = "StartracOut", definition = StartracOut.plot)


