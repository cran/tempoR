#' Write output
#'
#' Write the score data frame and plots of actual vs. predicted age for all reported gene sets to disk.
#'
#' @param results The results object returned by \code{\link{tempo.runInstance}}.
#' @param filename The prefix for the filename.  Writes filename.pdf and filename.table (a tab-delimited text file) for all reported gene sets
#' @examples
#' # In this example, text output will be written to tempfile.table and
#' # plots printed to tempfile.pdf
#' data("gse32472ExampleTempoResults")
#' tempo.writeOutput(gse32472ExampleTempoResults, tempfile(tmpdir = tempdir()))
#' @export
tempo.writeOutput = function(results, filename) {
    tempo.mkplots(results, sprintf("%s%s",filename,".pdf"))
    utils::write.table(results$scores[results$reported,],file=sprintf("%s%s",filename,".table"),quote=FALSE,sep="\t")
}

#' Make a plot for a specified gene sets
#'
#' Generate a plot of actual vs. predicted age, broken down by control sample or test sample, for a single specified gene set
#'
#' @param results The results object returned by \code{\link{tempo.runInstance}}.or \code{\link{tempo.run}}
#' @param geneset The name of a gene set
#' @examples
#' data("gse32472ExampleTempoResults")
#' tempo.mkplot(gse32472ExampleTempoResults, "cochlea development")
#' @export
tempo.mkplot = function(results,geneset){
  score = results$scores[geneset,"score"]
  ctrlResponse = results$Y[results$ctrl]
  testResponse = results$Y[results$test]
  ctrlPred = results$pred[geneset,results$ctrl]
  testPred = results$pred[geneset,results$test]
  yRange = c(min(ctrlResponse, testResponse)-1,max(c(ctrlResponse,testResponse)+1))
  graphics::par(mfrow=c(1,2),oma=c(4,4,2,0),mar=c(0,0,1,1))
  graphics::plot(ctrlResponse, ctrlPred, xlim = yRange, ylim = yRange, col = "red", main = "Control",axes=FALSE,asp=1)
  graphics::axis(side=1)
  graphics::axis(side=2)
  graphics::box()
  graphics::abline(0,1,lty=2)
  graphics::plot(testResponse, testPred, xlim = yRange, ylim = yRange, col = "blue", main = "Test",axes=FALSE,asp=1)
  graphics::axis(side=1)
  graphics::axis(side=2,labels=FALSE)
  graphics::box()
  graphics::abline(0,1,lty=2)
  graphics::title(main=sprintf("%s (Score: %f)", geneset, score), xlab = "Actual Age", ylab = "Predicted Age", outer=TRUE)
}

tempo.mkplots = function(results, fName) {
    grDevices::pdf(fName, width=9.75, height=6)
    for (geneset in results$reported) {
        tempo.mkplot(results,geneset)
    }
    grDevices::dev.off()
}
