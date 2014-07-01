\name{Polyfit-package}
\alias{Polyfit-package}
\alias{Polyfit}
\docType{package}
\title{
Polyfit add-on to DESeq and edgeR}
\description{
implementation the Polyfit add-on to DESeq and edgeR described in the paper 
"Improved error estimates for the analysis of differential expression from RNA-seq data"
}
\details{
\tabular{ll}{
Package: \tab Polyfit\cr
Type: \tab Package\cr
Version: \tab 0.99.2\cr
Date: \tab 2014-06-16\cr
License: \tab GPL(>=2)\cr
}
Polyfit is an add-on to the negative-binomial based packages \code{\link{edgeR}} and DESeq for two-class detection of differential expression which ensures the p-value distribution is uniform over the interval [0, 1] for data satisfying the null hypothesis of no differential expression.  The first component is the pair of functions \code{\link{pfExactTest}} and \code{\link{pfNbinomTest}} which replace the functions \code{\link{exactTest}} in edgeR and \code{\link{nbinomTest}} in DESeq respectively.  Their purpose is to smooth point singularities, particularly one at  \eqn{p = 1}, in the p-value distribution caused by calculating calculating p-values from a discrete distribution.  The output from either of these functions should then be passed to the second component, the function \code{link{levelPValues}}.  Its purpose is to apply a variant of the Storey-Tibshirani procedure to shift the p-values so that those corresponding to the null hypothesis have a unform distribution, and to calculate corresponding q-values (or 'adjusted p-values') for controlling errors via the false discovery rate.}
\author{
Conrad Burden

Maintainer:  conrad.burden@anu.edu.au
}
\references{
Burden, C.J., Qureshi, S. and Wilson, S.R. (2014). \emph{Error estimates for the 
analysis of differential expression from RNA-seq count data}, PeerJ PrePrints 2:e400v1. 

Robinson, M., McCarthy, D., and Smyth, G. (2010). \emph{edgeR: a Bioconductor package for differential expression analysis of digital gene expression data.} Bioinformatics, \bold{26}, 139-140.

Anders, S. and Huber, W. (2010). \emph{Differential expression analysis for sequence count data.} Genome Biology, \bold{11}(10), R106.
}
\keyword{ package }
\examples{
# Example using DESeq
cds <- makeExampleCountDataSet()
cds <- estimateSizeFactors( cds )
cds <- estimateDispersions( cds )
nbTPolyfit <- pfNbinomTest( cds, "A", "B" )

lP <- levelPValues(nbTPolyfit$pval)
pvalTab <- cbind(origPval=nbTPolyfit$pval, correctedPval=lP$pValueCorr, qval=lP$qValueCorr)
cat("\n Original and corrected P-values from DESeq \n")
head(pvalTab)

# Example using edgeR
countTable <- counts(cds)
nGenes <- dim(countTable)[1]
nCols <- dim(countTable)[2]
d <- DGEList(counts= countTable, group=c(1,1,2,2,2), lib.size=rep(nGenes, nCols))
 d <- calcNormFactors(d)
 d <- estimateCommonDisp(d)
 d <- estimateTagwiseDisp(d)
dePolyfit <- pfExactTest(d)

lP <- levelPValues(dePolyfit$table$PValue)
pvalTab <- cbind(origPval=dePolyfit$table$PValue, correctedPval=lP$pValueCorr, qval=lP$qValueCorr)
cat("\n Original and corrected P-values from edgeR \n")
head(pvalTab)
}
