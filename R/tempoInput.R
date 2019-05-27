#' Load a Gene Matrix Transposed formatted file
#'
#' \code{loadGMT} loads a Gene Matrix Transposed formatted file from a text file to the data structure used by TEMPO.  In the
#' .gmt format, each row represents a gene set, with tab delimited columns.  The first column is a gene set name,
#' the second columns and optional description, and the remaining columns contain gene ids for each gene in the
#' gene set.  The format is also described at
#' \href{http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats}{the BROAD site}.
#'
#' @param target a string indicating the location of the .gmt file
#' @return a list indexed by gene set name of lists of gene ids
#' @examples
#' # An example collection of gene sets is included in the package in .gmt format
#' exampleGeneSetsPath = file.path(path.package("tempoR"),"dflatExample.gmt")
#' exampleGeneSets = loadGMT(exampleGeneSetsPath)
#' @export
loadGMT = function(target) {
    fc = file(target)
    aList = strsplit(readLines(fc), "\t")
    close(fc)
    nms = sapply(aList,function(x){return(x[1])})
    out = sapply(aList,function(x){return(x[3:length(x)])})
    names(out) = nms
    return(out)
}

#' Load a Gene Cluster Text formatted file
#'
#' \code{loadGCT} loads a Gene Cluster Text formatted file from a text file to the data structure used by TEMPO. A .gct
#' file is organized as described at \href{http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats}{the BROAD site}.
#'
#' @param target a string indicating the location of the .gct file
#' @return a matrix with sample ids as row names and gene ids as column names
#' @examples
#' # An example gene expression data set is included in the package in .gct format
#' exampleDataPath = file.path(path.package("tempoR"),"gse32472Example.gct")
#' exampleData = loadGCT(exampleDataPath)
#' @export
loadGCT = function(target) {
    dat = utils::read.table(file=target, skip=2,header=TRUE,sep="\t")
    rownames(dat) = dat[,1]
    dat = dat[,3:dim(dat)[2]]
    return(t(as.matrix(dat)))
}

#' Load a categorical or continuous cls formatted file.
#'
#' \code{loadCLS} loads a categorical or continuous cls formated file.  Both file formats are described at \href{http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats}{the BROAD site}.
#'
#' @param target a string indicating the location of the .cls file. loadCLS will automatically determine whether the specified file is categorical or continuous.
#' @param sampleNames a list of sample names that correspond to each each entry in the .cls file.
#' @return a list indexed by sample id containing the specified categorical or continuous variable.
#' @examples
#' # The .cls format does not include sample names, so it is necesary to have some already.
#' # In this case, we get them from the .gct file
#' sampleNames = names(loadGCT(file.path(path.package("tempoR"),"gse32472Example.gct")))
#'
#' # Example continouous and categorical .cls files are included in the package
#' exampleAgesPath = file.path(path.package("tempoR"),"gse32472Example.age.cls")
#' exampleAges = loadCLS(exampleAgesPath, sampleNames)
#'
#' exampleClassesPath = file.path(path.package("tempoR"),"gse32472Example.phen.cls")
#' exampleClasses = loadCLS(exampleClassesPath, sampleNames)
#' @export
loadCLS = function(target, sampleNames) {
    fc = file(target)
    l = readLines(fc)
    close(fc)
    out = c()
    if(l[1] == "#numeric") {
        out = sapply(strsplit(l[3]," |\t")[[1]],as.numeric)
    } else {
        nms = strsplit(l[2]," ")[[1]]
        cls = strsplit(l[3]," |\t")[[1]]
        nminds = sapply((0:length(nms)-2),as.character)
        out = c()
        for(i in cls) {
            if(i %in% nminds){
              out = c(out,nms[as.integer(i)+2])
            } else {
                out = c(out,i)
            }
        }
    }
    names(out) = sampleNames
    return(out)
}
