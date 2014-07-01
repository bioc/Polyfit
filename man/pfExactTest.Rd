\name{pfExactTest}
\alias{pfExactTest}
\alias{pfExactTestDoubleTail}
\title{Polyfit extensions to the edgeR functions exactTest() and exactTestDoubleTail()
}

\description{Polyfit extension of the edgeR function \code{\link{exactTest}} and \code{\link{exactTestDoubleTail}} which perform exact genewise exact tests for differences in the means between two groups of negative-binomially distributed counts.}

\usage{
pfExactTest(object, pair=1:2, dispersion="auto", rejection.region="doubletail",
          big.count=900, prior.count.total=0.5)
pfExactTestDoubleTail(y1, y2, dispersion = 0, big.count = 900)
}

\arguments{
\item{object}{an object of class \code{\link[edgeR:DGEList-class]{DGEList}}.}

\item{pair}{vector of length two, either numeric or character, providing the pair of groups to be compared; if a character vector, then should be the names of two groups (e.g. two levels of \code{object$samples$group}); if numeric, then groups to be compared are chosen by finding the levels of \code{object$samples$group} corresponding to those numeric values and using those levels as the groups to be compared; if \code{NULL}, then first two levels of \code{object$samples$group} (a factor) are used. Note that the first group listed in the pair is the baseline for the comparison---so if the pair is \code{c("A","B")} then the comparison is \code{B - A}, so genes with positive log-fold change are up-regulated in group B compared with group A (and vice versa for genes with negative log-fold change).}

\item{dispersion}{either a numeric vector of dispersions or a character string indicating that dispersions should be taken from the data object.
If a numeric vector, then can be either of length one or of length equal to the number of tags.
Allowable character values are \code{"common"}, \code{"trended"}, \code{"tagwise"} or \code{"auto"}.
Default behavior (\code{"auto"} is to use most complex dispersions found in data object.}

\item{rejection.region}{type of rejection region for two-sided exact test.  Possible values are \code{"doubletail"}, \code{"smallp"} or \code{"deviance"}.}

\item{big.count}{count size above which asymptotic beta approximation will be used.}

\item{prior.count.total}{total prior count used to shrink log-fold-changes. Larger values produce more shrinkage.}

\item{y1}{numeric matrix of counts for the first the two experimental groups to be tested for differences.
Rows correspond to genes or transcripts and columns to libraries.
Libraries are assumed to be equal in size - e.g. adjusted pseudocounts from the output of \code{\link{equalizeLibSizes}}.}

\item{y2}{numeric matrix of counts for the second of the two experimental groups to be tested for differences.
Rows correspond to genes or transcripts and columns to libraries.
Libraries are assumed to be equal in size - e.g. adjusted pseudocounts from the output of \code{\link{equalizeLibSizes}}. Must have the same number of rows as \code{y1}.}
}

\value{
\code{exactTest} produces an object of class \code{DGEExact} containing the following components:
\item{table}{data frame containing columns for the log2-fold-change, \code{logFC}, the average log2-counts-per-million, \code{logCPM}, and the two-sided p-value \code{PValue}}
\item{comparison}{character vector giving the names of the two groups being compared} 
\item{genes}{optional data frame containing annotation for transcript; taken from \code{object}}

The low-level functions, \code{exactTestDoubleTail} etc, produce a numeric vector of genewise p-values, one for each row of \code{y1} and \code{y2}.
}

\details{
These functions have the same behaviour as the edgeR functions \code{\link{exactTest}} and \code{\link{exactTestDoubleTail}}, except that the `flagpole' in the P-value histogram, particularly at \eqn{p = 1} is redistributed using the function \code{\link{twoSidedPValueFromDiscrete}.}  
}

\references{
Burden, C.J., Qureshi, S. and Wilson, S.R. (2014). \emph{Error estimates for the 
analysis of differential expression from RNA-seq count data}, PeerJ PrePrints 2:e400v1. 

Robinson, M., McCarthy, D., and Smyth, G. (2010). \emph{edgeR: a Bioconductor package for differential expression analysis of digital gene expression data.} Bioinformatics, \bold{26}, 139-140.
}

\author{Conrad Burden, conrad.burden@anu.edu.au, based on existing software by Mark Robinson, Davis McCarthy, Gordon Smyth
}

\examples{
# generate raw counts from NB, create list object
y <- matrix(rnbinom(4000,size=1/0.2,mu=10),nrow=1000,ncol=4)
d <- DGEList(counts=y, group=c(1,1,2,2), lib.size=rep(1000,4))

de <- exactTest(d, dispersion=0.2)
topTags(de)

dePolyfit <- pfExactTest(d, dispersion=0.2)
topTags(dePolyfit)

oldpar <- par(mfrow=c(1,2))
hist(de$table$PValue,breaks=seq(0,1,by=0.01), 
  				xlab="P-value", main="edgeR")
hist(dePolyfit$table$PValue,breaks=seq(0,1,by=0.01), 
					xlab="P-value", main="polyfit-edgeR")
par(oldpar)
#
}

\seealso{
\code{\link{equalizeLibSizes}}, \code{\link{binomTest}}
}

