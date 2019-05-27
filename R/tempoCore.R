#' The main method.
#'
#' For each gene set in genesets, train a model on the control samples of Y~X.  Generates scores for each gene set, performs permutation testing, and (optionally) writes a table of results and plots
#'
#' @param X a matrix with sample ids as row names and gene ids as column names.
#' @param Y a list indexed by sample ids, containing numerical values.
#' @param genesets a list of lists.  Outer list is indexed by gene set name, inner list contains all gene ids in a given gene set
#' @param phen a list indexed by sample ids, containing phenotypes.  If ctrl and test are null, the phenotype of the first non-NA entry in the list is assumed to be the control phenotype; all others are test phenotypes
#' @param ctrl a list of sample ids.  The list of control samples to use in scoring.  If train is null, these are also the training samples.  Used only if phen is null.
#' @param test a list of sample ids.  The list of test samples to use in scoring.  Used only if phen is null.
#' @param train a list of sample ids.  The list of control samples to train models on.  If null, samples ids in ctrl are used. Used only if phen is null.
#' @param output optional.  A prefix to write table of output and plots for each reported gene set
#' @param numPerms number of permutations to do in permutation testing.  Defaults to 500
#' @param validation type of validation to do.  Defaults to leave one out, which generates deterministic results.  "CV" performs 10-fold cross-validation, which is significantly faster but generates non-deterministic results.
#' @param minGsSize minimum acceptable size for a gene set, considering only features which exist in \code{colnames(X)}
#' @param nCores number of thread to spawn for permutation testing.  This should likely be set to some number less than or equal to
#' the number of cores on your machine.  If nCores is less than 0, nCores will be set to the return value of \code{\link[parallel]{detectCores}}.
#' @param pCutoff report only gene sets with a p-value below this cutoff
#' @param fdrCutoff report only gene sets with a FDR below this cutoff
#' @param pMseCutoff report only gene sets with the p-value of the control mean squared error below this cutoff
#' @return the output from \code{\link{tempo.runInstance}} annotated with significance for each gene set
#' @examples
#' data("dflatExample")
#' data("gse32472Example")
#'
#' \donttest{
#' # This runs a simple TEMPO analysis on the example data set with default settings
#' # (with the exception of nCores, which will instead be automatically set to a suitable
#' # value) and saves the output in a two temporary files.
#' # Note that running this example may take several minutes.
#' results = tempo.run(phen=gse32472Example$bpd,
#'     genesets=dflatExample,
#'     X=gse32472Example$data,
#'     Y=gse32472Example$age,
#'     output=tempfile(tmpdir = tempdir()),
#'     nCores=-1)
#'
#' # If phen is used, the first item in the list is assumed to the control phenotype
#' # and all other phenotypes test. Specifiy ctrl and test exactly for more control.
#' # Note that running this example may take several minutes.
#' results = tempo.run(ctrl=gse32472Example$ctrl,
#'     test=gse32472Example$test,
#'     genesets=dflatExample,
#'     X=gse32472Example$data,
#'     Y=gse32472Example$age,
#'     nCores=-1)
#'
#'
#' # If training models on a held out set of data is desired, train can be specified seperately
#' # Note that running this example may take several minutes.
#' results2 = tempo.run(train=gse32472Example$ctrl[1:10],
#'     ctrl=gse32472Example$ctrl[11:20],
#'     test=gse32472Example$test,
#'     genesets=dflatExample,
#'     X=gse32472Example$data,
#'     Y=gse32472Example$age,
#'     nCores=-1)
#' }
#'
#' # Reporting thresholds, number of permutations, and number of CPU cores used can all be changed.
#' # This command is suitable for demonstration purposes, but significance values will not be
#' # meaningful.
#' results3 = tempo.run(phen=gse32472Example$bpd,
#'     genesets=dflatExample,X=gse32472Example$data,
#'     Y=gse32472Example$age,output=tempfile(tmpdir = tempdir()),
#'     numPerms=2,nCores=2,pCutoff=1,fdrCutoff=2,pMseCutoff = 1)
#' @export
tempo.run = function(X, Y, genesets, phen=NULL, ctrl=NULL, test=NULL, train=NULL, output = "", numPerms = 500, validation = "LOO", minGsSize=2, nCores=24,pCutoff=.05,fdrCutoff=.25,pMseCutoff=.05) {

    if(!is.null(phen)){
        phen = stats::na.omit(phen)
        ctrl = names(phen[phen == phen[[1]]])
        test = names(phen[phen != phen[[1]]])
        train = NULL
    }

    # remove any genes not in X from gene sets
    gs = c()
    for(i in names(genesets)){
        path = intersect(genesets[[i]],colnames(X))
        if(length(path)>=minGsSize){gs[[i]] = path}
    }
    genesets = gs

    # generate model scores for each gene set
    results = tempo.runInstance(X, Y, genesets, ctrl, test, train, validation=validation)

    # permutation testing to generate p and fdr for each gene set
    results = tempo.permutationTest(results, X, Y, genesets, ctrl, test, train, numPerms, nCores,pMseCutoff=pMseCutoff)
    results$reported = tempo.findReported(results$scores,pCutoff=pCutoff,fdrCutoff=fdrCutoff,pMseCutoff=pMseCutoff)

    if(output != ""){
        tempo.writeOutput(results,output)
    }
    return(results)
}

#' Build models for all pathways using the control data and test on the test population.
#'
#' @param X a matrix with sample ids as row names and gene ids as column names.
#' @param Y a list indexed by sample ids, containing numerical values.
#' @param genesets a list of lists.  Outer list is indexed by gene set name, inner list contains all gene ids in a given gene set
#' @param ctrl a list of sample ids.  The list of control samples to use in scoring.
#' @param test a list of sample ids.  The list of test samples to use in scoring.
#' @param  train a list of sample ids.  The list of control samples to train models on.  If null, train on ctrl.
#' @param comps maximum number of components to use in the plsr model
#' @param validation "CV" for 10-fold cross-validation, "LOO" for leave-one-out cross-validation. "CV" is nondeterministic and should not be used where exactly reproducible results are important
#' @return a list with the following entries
#' \itemize{
#'   \item \code{ctrl} a list of the sample names used from the control set used for scoring
#'   \item \code{test} a list of the sample names used from the test set used for scoring
#'   \item \code{train} a list of the sample names used from the control set to train models
#'   \item \code{scores} a data frame with gene set ids as row names and a "ctrlMSE" and "score" entry for each gene set with the MSE of control age predictions in cross-validation and the calulated score, respectively
#'   \item \code{pred} a matrix with gene set labels, where the age prediction for sample j from the models for gene set i are at i,j
#'   \item \code{Y} the continuous variable of interest that models are built with respect to, e.g. age
#'   \item \code{genesets} the gene sets used in the analysis
#' }
#' @examples
#' data("dflatExample")
#' data("gse32472Example")
#'
#' # It is possible to run just the model-building and scoring functions and skip
#' # cross-validation.  This is not recommended for general use, but may be useful
#' # in some cases
#' results = tempo.runInstance(ctrl=gse32472Example$ctrl,
#'     test=gse32472Example$test,
#'     genesets=dflatExample,
#'     X=gse32472Example$data,
#'     Y=gse32472Example$age)
#' @export
tempo.runInstance = function(X, Y, genesets, ctrl, test, train=NULL, comps = 10, validation="CV") {
    sortedNames = sort(unique(c(train,ctrl,test)))

    ctrl.X = X[ctrl,]
    ctrl.Y = Y[ctrl]
    train.X = ctrl.X
    train.Y = ctrl.Y
    test.X = X[test,]
    test.Y = Y[test]

    if(!is.null(train)){
        train.X = X[train,]
        train.Y = Y[train]
    }

    results = c()

    results$pred = matrix(0, nrow = length(genesets), ncol = length(sortedNames), dimnames=list(names(genesets),sortedNames))
    scores = data.frame(matrix(0, nrow = length(genesets), ncol = 2, dimnames=list(names(genesets),c("ctrlMSE","score"))))
    for(i in 1:length(genesets)) {
        gs = intersect(genesets[[i]],colnames(X))
        if(length(gs)<=1){ next }

        gsName = names(genesets[i])
        newF = train.X[, gs]

        print(sprintf("%s: %s", i, gsName))

        m = pls::plsr(train.Y ~ newF , validation = validation)

        # Choose number of components to use by minimizing mean squared error
        compMSE = c()
        for(j in  1:m$ncomp) {compMSE[[j]] = sum((m$validation$pred[,,j] - train.Y)**2)/length(train.Y)}

        ncomp = which.min(compMSE[1:min(comps,m$ncomp)])

        if(is.null(train)){
            ctrlPred = m$validation$pred[,,ncomp]
        } else {
            newF3 = ctrl.X[,gs]
            ctrlPred = stats::predict(m, newF3)[,,ncomp]
            results$pred[i,train] = m$validation$pred[,,ncomp]
        }
        results$pred[i,ctrl] = ctrlPred

        newF2 = test.X[, gs]
        testPred = stats::predict(m, newF2)[,,ncomp]
        results$pred[i,test] = testPred


        scores[gsName,"ctrlMSE"] = compMSE[[ncomp]]
	      scores[gsName,"score"] <- tempo.computeScore(ctrl, c(ctrlPred, testPred), Y)
    }

    results$train = names(train.Y)
    results$ctrl = names(ctrl.Y)
    results$test = names(test.Y)
    results$scores = scores
    results$Y = Y
    results$genesets = genesets

    return(results)
}

#' Permutation testing for TEMPO
#'
#' Does permutation testing for the results of an instance of tempo run on the true class assignments.  Called by \code{\link{tempo.run}}.
#'
#' @param mainResult the data structure returned by tempo.runInstance when run using the true class assignments
#' @param X a matrix with sample ids as row names and gene ids as column names
#' @param Y a list indexed by sample ids, containing numerical values, e.g. ages
#' @param genesets a list of lists.  Outer list is indexed by gene set name, inner list contains all gene ids in a given gene set
#' @param ctrl a list of sample ids.  The list of control samples to use in scoring.
#' @param test a list of sample ids.  The list of test samples to use in scoring.
#' @param train a list of sample ids.  The list of control samples to train models on.  If null, train on ctrl.
#' @param numPerms number of permutations to do in permutation testing.  Defaults to 500
#' @param nCores number of thread to spawn for permutation testing.  This should likely be set to some number less than or equal to
#' the number of cores on your machine.  If nCores is less than 0, nCores will be set to the return value of \code{\link[parallel]{detectCores}}.
#' @param pMseCutoff gene sets with a p-value associated with the control mean squared error below this cuttof are not evaluated further.
#' @return mainResult, with the score data frame annotated with p-values for ctrlMSE ("pMSE") and score ("p") and fdr for the score p-value ("BH") for each gene set.
#' Note that "BH" will be 2 for any gene set where pMSE > pMseCutoff
#' @examples
#' data("dflatExample")
#' data("gse32472Example")
#' data("gse32472ExampleTempoResults")
#'
#' \donttest{
#' # Cross-validation can be run seperately once the initial model-building and scoring is complete.
#' # Note that running this example may take several minutes.
#' results = tempo.runInstance(ctrl=gse32472Example$ctrl,
#'     test=gse32472Example$test,
#'     genesets=dflatExample,
#'     X=gse32472Example$data,
#'     Y=gse32472Example$age)
#' results = tempo.permutationTest(mainResult=results,
#'     ctrl=gse32472Example$ctrl,
#'     test=gse32472Example$test,
#'     genesets=dflatExample,
#'     X=gse32472Example$data,
#'     Y=gse32472Example$age,
#'     nCores=-1)
#' }
#'
#' # Reporting thresholds, number of permutations, and number of CPU cores used can all be changed.
#' # This command is suitable for demonstration purposes, but significance values will not be
#' # meaningful.
#' results2 = tempo.permutationTest(mainResult=gse32472ExampleTempoResults,
#'     ctrl=gse32472Example$ctrl,
#'     test=gse32472Example$test,
#'     genesets=dflatExample,
#'     X=gse32472Example$data,
#'     Y=gse32472Example$age,
#'     numPerms=2,nCores=2,pMseCutoff = 1)
#' @export
tempo.permutationTest = function(mainResult, X, Y, genesets, ctrl, test, train=NULL, numPerms = 500, nCores=24, pMseCutoff=.05){
    if(nCores < 0){
      nCores = parallel::detectCores()
    }
    print(sprintf("Beginning permutation testing for %s permutations on %s cores.",numPerms,nCores))
    usedCores = nCores
    cl = parallel::makeCluster(usedCores)
    doParallel::registerDoParallel(cl, cores = usedCores)

    ctrlSize = length(ctrl)
    n = c(ctrl,test)

    '%dopar%' <- foreach::"%dopar%"

    results= c()
    results = foreach::foreach(i = 1:numPerms,.errorhandling="remove", .inorder = FALSE, .packages = c("pls"),  .export = c("tempo.runInstance","tempo.computeScore")) %dopar% {
      newGenesets = c()
      genes = colnames(X)
      for(gs in names(genesets)){
        newGenesets[[gs]] = sample(genes,length(genesets[[gs]]))
      }
      pls.perm = tempo.runInstance(X , Y, newGenesets, ctrl, test, train=train)
        gc()
        result = pls.perm # This ends up stored in results
    }
    parallel::stopCluster(cl)

    numPerms = length(results)
    pScores = matrix(0, nrow = length(genesets), ncol = numPerms, dimnames=list(names(genesets),NULL))
    pCtrlMSE = matrix(0, nrow = length(genesets), ncol = numPerms, dimnames=list(names(genesets),NULL))
    pCtrl = list()
    for(i in 1:numPerms){
        pScores[,i] = results[[i]]$scores[,"score"]
        pCtrlMSE[,i] = results[[i]]$scores[,"ctrlMSE"]
        #pPreds[,,i] = results[[i]]$pred
        pCtrl[[i]] = results[[i]]$ctrl
    }
    allPermScores = c(pScores)
    allTrueScores = c(mainResult$scores[,"score"])
    scalar = length(allTrueScores)/length(allPermScores)
    for(i in names(genesets)){
        pathPerm = pScores[i,]
        thisScore = mainResult$scores[i,"score"]
        p = (1+length(pathPerm[pathPerm >= thisScore]))/(numPerms+1)
        msePerm = pCtrlMSE[i,]
        mseTrue = mainResult$scores[i,"ctrlMSE"]
        mainResult$scores[i,"pMSE"] = (1+length(msePerm[msePerm <= mseTrue]))/(numPerms+1)
        mainResult$scores[i,"p"] = p
    }
    i1 = mainResult$scores[,"pMSE"] <= pMseCutoff
    mainResult$scores[,"BH"] = rep(2,length(rownames(mainResult$scores)))
    mainResult$scores[i1,"BH"] = stats::p.adjust(mainResult$scores[i1,"p"],method="BH")
    return(mainResult)
}

tempo.findReported = function(scoreDF,pCutoff=.05,fdrCutoff=.25,pMseCutoff=.05) {
    i2 = scoreDF["p"] <= pCutoff
    i3 = scoreDF["BH"] <= fdrCutoff
    i1 = scoreDF["pMSE"] <= pMseCutoff

    d = scoreDF[i1 & i2 & i3,]

    # return a list of names sorted by score
    sigNames = rownames(d)
    d = d[["score"]]
    names(d) = sigNames
    return(names(d[sort.list(d,decreasing=TRUE)]))
}

tempo.computeScore = function(ctrl, preds, Y, filterType = 0) {
    Y = Y[names(preds)]
    if (filterType==1){
        preds2 = preds[preds > 2*min(Y) - mean(Y) & preds < 2*max(Y) - mean(Y)]
    } else {
        preds2 = preds
    }
    nm = names(preds2)
    ctrl2 = intersect(nm, ctrl)
    test2 = setdiff(nm,ctrl2)

    ctrlErr = Y[ctrl2]-preds2[ctrl2]
    testErr = Y[test2]-preds2[test2]

    # Learn the MLE normal distribution on the control data
    ctrlMean = mean(ctrlErr)
    ctrlSD = stats::sd(ctrlErr)

    s1 = 0
    s2 = 0

	minval <- 10^-200

    # Calculate mean likelihood of seeing the control errors on the learned error model
    for(j in ctrlErr) {
        if(j > ctrlMean){
            p = stats::pnorm(j, mean=ctrlMean, sd=ctrlSD, lower.tail=F)
        } else {
            p = stats::pnorm(j, mean=ctrlMean, sd=ctrlSD)
        }
		s1 = s1 - log(max(p,minval)) #set a minimum so pnorm returning 0 doesn't carry through to produce infinite scores
    }
    s1 = s1/length(ctrlErr)

    # Calculate mean likelihood of seeing the test errors on the learned error model
    for(j in testErr) {
        if(j > ctrlMean){
            p = stats::pnorm(j, mean=ctrlMean, sd=ctrlSD, lower.tail=F)
        } else {
            p = stats::pnorm(j, mean=ctrlMean, sd=ctrlSD)
        }
		s2 = s2 - log(max(p,minval))
    }
    s2 = s2/length(testErr)

    return(s2/s1)
}
