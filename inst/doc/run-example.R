## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----load,eval=TRUE------------------------------------------------------
library("tempoR")
data("dflatExample")
data("gse32472Example")

## ----create, eval=FALSE--------------------------------------------------
#  dflatExample = loadGMT(paste0(path.package("tempoR"),"/dflatExample.gmt"))
#  gse32472Example = list()
#  gse32472Example$data = loadGCT(paste0(path.package("tempoR"),"/gse32472Example.gct"))
#  gse32472Example$age = loadCLS(paste0(path.package("tempoR"),"/gse32472Example.age.cls"),
#                                sampleNames=rownames(gse32472Example$data))
#  gse32472Example$bpd = loadCLS(paste0(path.package("tempoR"),"/gse32472Example.phen.cls"),
#                                sampleNames=rownames(gse32472Example$data))
#  
#  gse32472Example$ctrl = names(which(gse32472Example$bpd=="no"))
#  gse32472Example$test = names(which(gse32472Example$bpd!="no"))

## ----run0, eval=FALSE----------------------------------------------------
#  results = tempo.run(phen=gse32472Example$bpd,
#                      genesets=dflatExample,
#                      X=gse32472Example$data,
#                      Y=gse32472Example$age,
#                      numPerms=4,
#                      nCores=2,
#                      output=tempfile(tmpdir = tempdir()),
#                      pCutoff=1,
#                      fdrCutoff=2,
#                      pMseCutoff = 1)

## ----run2, eval=FALSE----------------------------------------------------
#  results = tempo.run(ctrl=gse32472Example$ctrl,
#                      test=gse32472Example$test,
#                      genesets=dflatExample,
#                      X=gse32472Example$data,
#                      Y=gse32472Example$age,
#                      numPerms=4,
#                      nCores=42
#                      output=tempfile(tmpdir = tempdir()),
#                      pCutoff=1,
#                      fdrCutoff=2,
#                      pMseCutoff = 1)

## ----runActual, message=FALSE, warning=FALSE, include=FALSE--------------
library("tempoR")
results = tempo.run(ctrl=gse32472Example$ctrl,test=gse32472Example$test,genesets=dflatExample,X=gse32472Example$data,Y=gse32472Example$age,numPerms=4,nCores=2,pCutoff=1,fdrCutoff=2,pMseCutoff = 1)
knitr::kable(results$scores[results$reported,])

## ----runTable, echo=FALSE------------------------------------------------
knitr::kable(results$scores[results$reported,])

## ----runPlot, fig.width=7, fig.height=4.5--------------------------------
tempo.mkplot(results,"cochlea development")

## ----runPlot2, fig.width=7, fig.height=4.5-------------------------------
tempo.mkplot(results,"regulation of bone remodeling")

## ----run3, eval=FALSE----------------------------------------------------
#  results2 = tempo.run(train=gse32472Example$ctrl[1:10],
#                       ctrl=gse32472Example$ctrl[11:20],
#                       test=gse32472Example$test,
#                       genesets=dflatExample,
#                       X=gse32472Example$data,
#                       Y=gse32472Example$age,
#                       numPerms=4,
#                       nCores=4)

