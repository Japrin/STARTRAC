
#' dispaly message with time stamp
#' @param msg characters; message to display
#' @export
loginfo <- function(msg) {
  timestamp <- sprintf("%s", Sys.time())
  msg <- paste0("[",timestamp, "] ", msg,"\n")
  cat(msg)
}

#' entropy of each row of the input matrix
#' @param x matrix;
mrow.entropy <- function(x)
{
  freqs <- sweep(x,1,rowSums(x),"/")
  H = - rowSums(ifelse(freqs>0,freqs* log2(freqs),0))
  return(H)
}

#' entropy of each column of the input matrix
#' @param x matrix;
mcol.entropy <- function(x)
{
  freqs <- sweep(x,2,colSums(x),"/")
  H = - colSums(ifelse(freqs>0,freqs* log2(freqs),0))
  return(H)
}

#' gini_simpson of each column of the input matrix
#' @param x matrix;
mcol.gini_simpson <- function(x)
{
  freqs <- sweep(x,2,colSums(x),"/")
  GS <- 1-colSums(freqs^2)
  return(GS)
}

#' Gini of the input vector
#' @param y vector;
Gini <- function(y)
{
  ### https://en.wikipedia.org/wiki/Gini_coefficient
  y <- y[y>0]
  y <- sort(y)
  n <- length(y)
  G <- 1/n * (n+1 - 2* sum( (n + 1 - 1:n)*y) / sum(y))
  return(G)
}


#' warpper function for Startrac analysis
#' @importFrom data.table dcast
#' @importFrom plyr ldply adply llply
#' @importFrom parallel makeCluster stopCluster
#' @importFrom RhpcBLASctl omp_set_num_threads
#' @importFrom doParallel registerDoParallel
#' @importFrom methods new slot
#' @importFrom methods slot<-
#' @param cell.data data.frame. Each line for a cell, and these columns as required: `Cell_Name`, `clone.id`, `patient`, `majorCluster`, `loc`
#' @param proj character. String used to annotate the project.
#' @param cores integer. number of core to be used. default: NULL.
#' @param n.perm integer. number of permutation will be performed. If NULL, no permutation. (default: NULL)
#' @param verbose integer. verbose level indicate how to include intermediate result (some Startrac objects)
#' @details run the Startrac pipeline
#' @return an list contains data.frame elements "cluster.data","pIndex.migr" and "pIndex.tran"
#' @export
#' @examples 
#' library("Startrac")
#' dat.file <- system.file("extdata/example.cloneDat.Zhang2018.txt",package = "Startrac")
#' in.dat <- read.table(dat.file,stringsAsFactors = FALSE,head=TRUE)
#' out <- Startrac.run(in.dat, proj="CRC", cores=2,verbose=FALSE)
#' 
Startrac.run <- function(cell.data, proj="CRC", cores=NULL,n.perm=NULL,verbose=0)
{
  ##tic("obj.proj")
  RhpcBLASctl::omp_set_num_threads(1)
  loginfo("initialize Startrac ...")
  obj.proj <- new("Startrac",cell.data,aid=proj,n.perm=n.perm,cores=cores)
  loginfo("calculate startrac index ...")
  obj.proj <- calIndex(obj.proj,cores=cores,n.perm=n.perm)
  loginfo("calculate pairwise index ...")
  obj.proj <- pIndex(obj.proj,cores=cores,n.perm=n.perm)
  if(!is.null(n.perm)){ 
    loginfo("get the significance")
    obj.proj <- getSig(obj.proj,obj.proj@cell.perm.data) 
  }else{
    obj.proj <- getSig(obj.proj,NULL) 
  }
  ##toc()
  
  obj.list <- NULL
  if(length(obj.proj@patient.size)>1)
  {
    loginfo("calculate indices of each patient ...")
    patient.vec <- names(obj.proj@patient.size[obj.proj@patient.size > 30])
    #cl <- makeCluster(if(is.null(cores)) (detectCores()-2) else cores)
    #registerDoParallel(cl)
    withCallingHandlers({
      obj.list <- llply(patient.vec,function(pid,cell.data){
        require("Startrac")
        obj <- new("Startrac",subset(cell.data,patient==pid),aid=pid)
        obj <- calIndex(obj)
        obj <- pIndex(obj,cores=1)
        obj <- getSig(obj,NULL)
        obj
      },cell.data=cell.data,.progress = "none",.parallel=F)
      names(obj.list) <- patient.vec
    },warning=function(w) {
      if(grepl("... may be used in an incorrect context:",conditionMessage(w)))
        ### strange bug, see https://github.com/hadley/plyr/issues/203
        invokeRestart("muffleWarning")
    })
    #stopCluster(cl)
  }
  loginfo("collect result")
  ret <- new("StartracOut",proj=proj)
  ## cluster index
  ret.slot.names <- c("cluster.data","pIndex.migr","pIndex.tran",
                      "cluster.sig.data","pIndex.sig.migr","pIndex.sig.tran")
#  if(!is.null(n.perm)){ 
#    ret.slot.names <- c(ret.slot.names,
#                        c("cluster.sig.data","pIndex.sig.migr","pIndex.sig.tran")) 
#  }
  for(v in ret.slot.names)
  {
    .tmp.list <- c(obj.proj,obj.list)
    names(.tmp.list)[1] <- proj
    slot(ret,v) <- ldply(.tmp.list,function(obj){ slot(obj,v) },.id=NULL)
    .tmp.list <- NULL
#    slot(ret, v) <- slot(obj.proj,v)
#    if(!is.null(obj.list)){
#      slot(ret, v) <- rbind(slot(ret, v),ldply(obj.list,function(obj){
#        slot(obj,v)
#      }))
#    }
  }
  if(verbose==1){
      ret@objects <- c(obj.proj)
  }else if(verbose==2){
      ret@objects <- c(obj.proj,obj.list)
  }
  loginfo("return")
  return(ret)
}

#' calculate Startrac.dist (tissue distribution preference)
#' @import data.table
#' @importFrom plyr aaply
#' @importFrom stats chisq.test
#' @param dat.tb data.frame. Each line for a cell, and these columns as required: `majorCluster`, `loc`
#' @param byPatient logical. whether calculate the index for each patient. (default: FALSE)
#' @param colname.cluster character. which column specify the cluster (default: "majorCluster")
#' @param colname.patient character. which column specify the patient  (default: "patient")
#' @param colname.tissue character. which column specify the tissue  (default: "loc")
#' @param method character. method to use, one of "chisq", "fisher", and "freq"  (default: "chisq")
#' @param min.rowSum integer. rows with rowSum <= this value will be filtered out (default: 0)
#' @details calculate Startrac.dist (tissue distribution preference).
#' @return an array full of R_{o/e} (method="chisq") or list with components of OR, p.value etc. from fisher.test (method="fisher")
#' @export
calTissueDist <- function(dat.tb,byPatient=F,colname.cluster="majorCluster",
							  colname.patient="patient",colname.tissue="loc",
                              method="chisq",min.rowSum=0)
{

    .table.fisher <- function(count.dist)
    {
        sum.col <- colSums(count.dist)
        sum.row <- rowSums(count.dist)
        count.dist.tb <- as.data.frame(unclass(count.dist))
        setDT(count.dist.tb,keep.rownames=T)
        count.dist.melt.tb <- melt(count.dist.tb,id.vars="rn")
        colnames(count.dist.melt.tb) <- c("rid","cid","count")
        count.dist.melt.ext.tb <- as.data.table(ldply(seq_len(nrow(count.dist.melt.tb)), function(i){
                               this.row <- count.dist.melt.tb$rid[i]
                               this.col <- count.dist.melt.tb$cid[i]
                               this.c <- count.dist.melt.tb$count[i]
                               other.col.c <- sum.col[this.col]-this.c
                               this.m <- matrix(c(this.c,
                                          sum.row[this.row]-this.c,
                                          other.col.c,
                                          sum(sum.col)-sum.row[this.row]-other.col.c),
                                        ncol=2)
                               res.test <- fisher.test(this.m)
                               data.frame(rid=this.row,
                                      cid=this.col,
                                      p.value=res.test$p.value,
                                      OR=res.test$estimate)
                               }))
        count.dist.melt.ext.tb <- merge(count.dist.melt.tb,count.dist.melt.ext.tb,
                                        by=c("rid","cid"))
        count.dist.melt.ext.tb[,p.adj:=p.adjust(p.value,"BH")]
        return(count.dist.melt.ext.tb)
    }

    if(method=="freq"){
        ncount.sampleID <- dat.tb[,.(N.tot=.N),by=c("sampleID","loc")]
        ncount.sampleID_mcls <- dat.tb[,.(N.mcls=.N),by=c("sampleID","loc","majorCluster")]
        ncount.sampleID_mcls <- melt(dcast(ncount.sampleID_mcls,sampleID+loc~majorCluster,
                                           value.var="N.mcls",fill=0),
                                     id.vars=c("sampleID","loc"),
                                     variable.name="majorCluster",
                                     value.name="N.mcls")
        ncount.sampleID_mcls <- merge(ncount.sampleID_mcls,ncount.sampleID,by=c("sampleID","loc"))
        ncount.sampleID_mcls[,freq.mcls:=N.mcls/N.tot]
        res <- ncount.sampleID_mcls[,{
                                        loc.vec <- unique(.SD$loc)
                                        o.tb <- as.data.table(ldply(loc.vec,function(xx){
                                            freq.x <- .SD[loc==xx,][["freq.mcls"]]
                                            #freq.y <- .SD[["freq.mcls"]]
                                            freq.y <- .SD[loc!=xx,][["freq.mcls"]]
                                            if(length(freq.x) >= 3)
                                            {
                                                oo.t <- t.test(freq.x,freq.y)
                                                oo.w <- wilcox.test(freq.x,freq.y)
                                                data.table(loc=xx,
                                                       logFC=log2(mean(freq.x)/mean(freq.y)),
                                                       p.value.t=oo.t$p.value,
                                                       p.value.w=oo.w$p.value)
                                            }else{
                                                NULL
                                            }
                                        }))
                                    },by=c("majorCluster")][order(majorCluster,loc),]
        res[,FDR.t:=p.adjust(p.value.t,"BH")]
        res[,FDR.w:=p.adjust(p.value.w,"BH")]

        res[,char.sig:=""]
        res[FDR.t < 0.01 & p.value.t < 0.05,char.sig:="\U2020"]
        res[FDR.t < 0.05,char.sig:="\U2731"]
        res[FDR.t < 0.01,char.sig:="\U2731\U2731"]
        #res[FDR.t < 0.05,char.sig:="*"]
        #res[FDR.t < 0.01,char.sig:="**"]
        dist.charSig.tb <- dcast(res,majorCluster~loc,value.var="char.sig")
        dist.logFC.tb <- dcast(res,majorCluster~loc,value.var="logFC")
        dist.FDR.tb <- dcast(res,majorCluster~loc,value.var="FDR.t")
        ret <- list("res"=res,
                    "dist.logFC.tb"=dist.logFC.tb,
                    "dist.FDR.tb"=dist.FDR.tb,
                    "dist.charSig.tb"=dist.charSig.tb)


    }else{
        if(byPatient==F){
            N.o <- table(dat.tb[[colname.cluster]],dat.tb[[colname.tissue]])
            if(method=="chisq"){
                res.chisq <- chisq.test(N.o)
                R.oe <- (res.chisq$observed)/(res.chisq$expected)
                ret <- R.oe
            }else if(method=="fisher"){
                count.dist <- N.o[rowSums(N.o) > min.rowSum,,drop=F]
                count.dist.melt.ext.tb <- .table.fisher(count.dist)
                dist.p.tb <- dcast(count.dist.melt.ext.tb,rid~cid,value.var="p.value")
                dist.padj.tb <- dcast(count.dist.melt.ext.tb,rid~cid,value.var="p.adj")
                dist.OR.tb <- dcast(count.dist.melt.ext.tb,rid~cid,value.var="OR")
                ret <- list("count.dist"=count.dist.melt.ext.tb,
                            "p.tb"=dist.p.tb,
                            "padj.tb"=dist.padj.tb,
                            "OR.tb"=dist.OR.tb)
            }
        }else{
            N.o.byPatient <- table(dat.tb[[colname.patient]],
                                   dat.tb[[colname.cluster]], dat.tb[[colname.tissue]])
            ret <- aaply(N.o.byPatient,1,function(x){
                             if(method=="chisq"){
                                res.chisq <- chisq.test(x)
                                return((res.chisq$observed)/(res.chisq$expected))
                             }else{
                                res.fisher <- .table.fisher(x)
                                return(dcast(res.fisher,rid~cid,value.var="OR"))
                             }
                    })
        }
    }
	return(ret)
}

#' plot Startrac.dist (tissue distribution preference)
#' @import data.table
#' @importFrom sscVis plotMatrix.simple
#' @importFrom grid gpar
#' @param OR.mtx matrix. OR data
#' @param k integer. for row-clustering. (default: 2)
#' @param method.distance character. for row-clustering. (default: "cosine")
#' @param do.hclust logical. for row-clustering. (default: TRUE)
#' @param out.prefix character. out.prefix  (default: NULL)
#' @param OR.max double. maximum of OR. ORs > this value will be set to this value (default: 3)
#' @param OR.min double. minimum of OR. ORs < this value will be set to this value (default: 0)
#' @param col.rid character. column indicating row ID (default: "rid")
#' @param col.ht vector;
#' @param exp.name character. legend title (default: expression(italic(OR))
#' @param p.tb data.table. p.value table. A column indicated by col.rid is required. (default: NULL)
#' @param charSig.tb data.table. charSig table. A column indicated by col.rid is required. (default: NULL)
#' @param pdf.width double. pdf width  (default: 5.5)
#' @param pdf.height double. pdf height  (default: 10)
#' @param ... parameters passed to sscVis:::plotMatrix.simple().
#' @details plot Startrac.dist (tissue distribution preference).
#' @return object returned by sscVis:::plotMatrix.simple()
#' @export
plotTissueDist <- function(OR.mtx,k=2,method.distance="cosine",do.hclust=T,
                            out.prefix=NULL,
                            OR.max=3,
                            OR.min=0,
                            col.rid="rid",
                            col.ht=circlize::colorRamp2(c(0, 1, 3), viridis::viridis(3)),
                            exp.name=expression(italic(OR)),
                            p.tb=NULL,
                            charSig.tb=NULL,
                            pdf.width = 5.5, pdf.height = 10,...)
{
    ### show.number
    #OR.mtx.tmp.txt <- apply(OR.mtx,2,function(x){ ifelse(x > OR.max, sprintf(">%1.0f",OR.max),sprintf("%1.2f",x)) })
    #OR.mtx.tmp.txt <- apply(OR.mtx.tmp.txt,2,function(x){ ifelse(x < OR.min, sprintf("<%1.0f",OR.min),sprintf("%1.2f",x)) })
    OR.mtx.tmp.txt <- apply(OR.mtx,2,function(x){ ifelse(x > OR.max,
                                                         sprintf(">%1.0f",OR.max),
                                                         ifelse(x < OR.min,sprintf("<%1.0f",OR.min),sprintf("%1.2f",x))
                                                         ) })
    rownames(OR.mtx.tmp.txt) <- rownames(OR.mtx)
    if(!is.null(p.tb)){
        p.mtx <- as.matrix(p.tb[,-c(col.rid),with=F])
        rownames(p.mtx) <- p.tb[[col.rid]]
        p.mtx <- p.mtx[rownames(OR.mtx.tmp.txt),colnames(OR.mtx.tmp.txt),drop=F]
        #p.mtx.txt <- apply(p.mtx,2,function(x){ ifelse(x < 0.01,"*","") })
        p.mtx.txt <- apply(p.mtx,2,function(x){ ifelse(x < 0.05,ifelse(x < 0.01,"**","*"),"") })
        p.mtx.txt[is.na(p.mtx.txt)] <- ""
        show.tmp.txt <- matrix(paste0(OR.mtx.tmp.txt,p.mtx.txt),nrow=nrow(OR.mtx.tmp.txt))
        rownames(show.tmp.txt) <- rownames(OR.mtx.tmp.txt)
        colnames(show.tmp.txt) <- colnames(OR.mtx.tmp.txt)
    }else if(!is.null(charSig.tb)){
        charSig.mtx <- as.matrix(charSig.tb[,-c(col.rid),with=F])
        rownames(charSig.mtx) <- charSig.tb[[col.rid]]
        charSig.mtx <- charSig.mtx[rownames(OR.mtx),colnames(OR.mtx),drop=F]
        show.tmp.txt <- charSig.mtx
    }else{
        show.tmp.txt <- OR.mtx.tmp.txt
    }

    ### clamp 
    OR.mtx.tmp <- OR.mtx
    OR.mtx.tmp[OR.mtx.tmp > OR.max] <- OR.max
    OR.mtx.tmp[OR.mtx.tmp < OR.min] <- OR.min
    OR.hclust.row <- run.cutree(OR.mtx.tmp,k=k,method.distance=method.distance,method.hclust="ward.D2")
    OR.hclust.row$branch <- dendextend::set(OR.hclust.row$branch,"branches_lwd", 2)
    #### ComplexHeatmap_2.2.0 work fine!
    #### ComplexHeatmap 2.7.1.1011  weird behaviour
    sscVis:::plotMatrix.simple(OR.mtx.tmp,
                             col.ht=col.ht,
                             out.prefix=sprintf("%s.tissue.dist.rClust.withDend",out.prefix),
                             returnHT = T,
                             show.number=show.tmp.txt,
                             show.dendrogram=do.hclust,
                             clust.row=if(do.hclust) OR.hclust.row$branch else FALSE,
                             row_dend_width = unit(1.5, "cm"),
                             #par.legend=list(color_bar = "discrete",at=seq(0,4,0.5)),
                             #waterfall.row=T,par.warterfall = list(score.alpha = 2,do.norm=T),
                             exp.name=exp.name,
                             z.hi=OR.max,
                             z.lo=OR.min,
                             mytitle = "Tissue Distribution",
                             #palatte=rev(brewer.pal(n = 7,name = "RdYlBu")),
                             #palatte=viridis::viridis(7),
                             par.heatmap=list(cex.row=1.5,
                                              row_names_gp=grid::gpar(
                                                #col=col.row.man,
                                                fontsize=10)),
                             pdf.width = pdf.width, pdf.height = pdf.height,...)
  
}


#' calculate the log likelihood ratio of each clonotype, accounting the uncertainty of cluster label assignment
#' @import data.table
#' @importFrom plyr llply
#' @param cellData data.table. Each line for a cell, and these columns as required: `majorCluster`
#' @param prob.mat matrix Lines for cells, and columns for assignment probability of the majorClusters. row names and column names are required 
#' @param cloneID.x character. If specified, it will subset the cellData using clonotype id  (default: NULL)
#' @param verbose logical. control the returned data types. (default: FALSE)
#' @return a list (verbose==T) or data.table (verbose==F)
#' @export
calCloneLLR <- function(cellData,prob.mat,cloneID.x=NULL,verbose=F)
{
    if(!is.null(cloneID.x)){
	cellData <- cellData[cloneID==cloneID.x,]
    }

    cellData[,majorCluster:=as.character(majorCluster)]
    cellData <- cellData[order(majorCluster),]
    cellData.N.tb <- cellData[,.N,by="majorCluster"][order(majorCluster),]
####	if(!is.null(mcls.range)){
####	    patch.N.tb <- data.table(majorCluster=setdiff(mcls.range,cellData.N.tb$majorCluster), N=0)
####	    cellData.N.tb <- rbind(cellData.N.tb,patch.N.tb)
####	}
    cellData.N.vec <- structure(cellData.N.tb$N,names=cellData.N.tb$majorCluster)
    cellData.f.vec <- cellData.N.vec/sum(cellData.N.vec)
    G.possible.f <- do.call(c,llply(seq_len(length(cellData.f.vec)),function(l){
	combn(cellData.f.vec,l,simplify=F)
    }))
    names(G.possible.f) <- sapply(G.possible.f,function(x){ paste0(names(x),collapse=":")  })

    LL <- do.call(c,llply(names(G.possible.f),function(g){
	g.f <- G.possible.f[[g]]
	## L == P(D|G==g)
	## LL == log(P(D|G==g)) ==sum(log(P(dk|G==g)))
	LL.g <- sum(apply(cellData,1,function(dk){
			    prob.k <- prob.mat[ dk[["Cell_Name"]], ]
			    log10(g.f %*% prob.k[names(g.f)])
			}))
	return(LL.g)
    },.parallel=F))
    names(LL) <- names(G.possible.f)
    LL <- sort(LL,decreasing=T)
    LLR <- unname(LL[1]-LL[2])
    G.best <- names(LL)[1]
    G.obs <- paste0(names(cellData.f.vec),collapse=":")
    if(verbose){
	return(list("LL"=LL,"LLR"=LLR,"G.best"=G.best,"G.obs"=G.obs))
    }else{
	return(data.table("cloneID"=cloneID.x,"LLR"=LLR,"G.best"=G.best,"G.obs"=G.obs))
    }
}


