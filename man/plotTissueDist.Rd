% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{plotTissueDist}
\alias{plotTissueDist}
\title{plot Startrac.dist (tissue distribution preference)}
\usage{
plotTissueDist(
  OR.mtx,
  k = 2,
  method.distance = "cosine",
  do.hclust = T,
  out.prefix = NULL,
  OR.max = 3,
  OR.min = 0,
  col.rid = "rid",
  col.ht = circlize::colorRamp2(c(0, 1, 3), viridis::viridis(3)),
  exp.name = expression(italic(OR)),
  p.tb = NULL,
  charSig.tb = NULL,
  mytitle = "Tissue Distribution",
  pdf.width = 5.5,
  pdf.height = 10,
  ...
)
}
\arguments{
\item{OR.mtx}{matrix. OR data}

\item{k}{integer. for row-clustering. (default: 2)}

\item{method.distance}{character. for row-clustering. (default: "cosine")}

\item{do.hclust}{logical. for row-clustering. (default: TRUE)}

\item{out.prefix}{character. out.prefix  (default: NULL)}

\item{OR.max}{double. maximum of OR. ORs > this value will be set to this value (default: 3)}

\item{OR.min}{double. minimum of OR. ORs < this value will be set to this value (default: 0)}

\item{col.rid}{character. column indicating row ID (default: "rid")}

\item{col.ht}{vector;}

\item{exp.name}{character. legend title (default: expression(italic(OR))}

\item{p.tb}{data.table. p.value table. A column indicated by col.rid is required. (default: NULL)}

\item{charSig.tb}{data.table. charSig table. A column indicated by col.rid is required. (default: NULL)}

\item{pdf.width}{double. pdf width  (default: 5.5)}

\item{pdf.height}{double. pdf height  (default: 10)}

\item{...}{parameters passed to sscVis:::plotMatrix.simple().}
}
\value{
object returned by sscVis:::plotMatrix.simple()
}
\description{
plot Startrac.dist (tissue distribution preference)
}
\details{
plot Startrac.dist (tissue distribution preference).
}
