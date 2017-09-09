##' @importFrom biovizBase flatGrl
##' @importFrom GenomicAlignments readGAlignments GAlignments GAlignmentPairs
##' @importFrom GenomicRanges GRangesList
##' @import  methods
##' @import  GenomeInfoDb
##' @import dplyr
##' @import magrittr
NULL

##' @rdname ChIPdata
##' @export
setClass("ChIPdata",
         representation = representation(
           readsSE = "GAlignments",
           readsPE = "GAlignmentPairs",
           fwdCover = "RleList",
           revCover = "RleList",
           isPE = "logical",
           metadata = "list"
         ),
         prototype = prototype(
           readsSE = GAlignments(),
           readsPE = GAlignmentPairs(first = GAlignments(),last = GAlignments()),
           isPE = FALSE,
           metadata = list()
         ))


##' ChIPdata object and constructors
##' 
##' \code{ChIPdata} is used to contain either SE or PE reads, and the idea
##' is to contain some of the most commonly used functions in a 
##' ChIP-seq/exo/nexus analysis.
##' 
##' @param bamfile a string indicating the location of the aligned reads bamfile
##' 
##' @param reads either a \code{GAlignments} or a \code{GAlignmentPairs} object with
##' the aligned reads of the ChIP experiment
##' 
##' @param isPE a logical value indicating if the reads are Single-End (SE) or Paired-End (PE). 
##' In the case that both \code{reads} and \code{isPE} are both used, then the \code{isPE} value
##' will be determined from \code{reads} and raise a warning
##' 
##' @return a \code{ChIPdata} object
##' 
##' @aliases ChIPdata, ChIPdata-class
##' @docType class
##' 
##' @examples
##' file <- system.file("extdata","example","encode_K562_Ctcf_first3chr_Rep1.sort.bam",package = "ChIPUtils")
##' ChIPdata(bamfile = file,isPE = FALSE)
##' 
##' reads <- readGAlignments(file,param = NULL)
##' ChIPdata(reads = reads)
ChIPdata <- function(bamfile = NULL,reads = NULL,isPE = NULL)
{
  if(is.null(bamfile) & is.null(reads)){
    stop("Both bamfile and reads are empty")
  }
  
  if(!is.null(bamfile) & !is.null(reads)){
    stop("Both bamfile and reads are not-empty, can't choose which one to use")
  }
  
  emptySE <- GAlignments()
  emptyPE <- GAlignmentPairs(first = emptySE,last = emptySE)
  
  if(!is.null(bamfile)){
    stopifnot(file.exists(bamfile),
              is.logical(isPE))  
    
    if(!isPE){
      reads <- readGAlignments(bamfile,param = NULL)
      sl <- seqlengths(reads)
      gl <- as(reads,"GRangesList")
      gl <- resize(gl,width = 1)
    }else{
      reads <- readGAlignmentPairs(bamfile,param = NULL)
      sl <- seqlengths(reads)
      gl <- as(reads,"GRangesList")
    }
    
  }else{
    possibleClass <- c(class(emptySE),class(emptyPE))
    stopifnot(class(reads) %in% possibleClass)

    isPE <- class(reads) == class(emptyPE)
    
    if(!isPE){
      sl <- seqlengths(reads)
      gl <- as(reads,"GRangesList")
      gl <- resize(gl,width = 1)
    }else{
      sl <- seqlengths(reads)
      gl <- as(reads,"GRangesList")
    }
  
  }

  fwdCover <- coverage(
    BiocGenerics::subset(gl,strand(gl) == "+"),
    width = sl)
  revCover <- coverage(
    BiocGenerics::subset(gl,strand(gl) == "-"),
    width = sl)
  
  if(!isPE){
    chipdata <- new("ChIPdata",
                    readsSE = reads,
                    readsPE = emptyPE,
                    isPE = FALSE,
                    fwdCover = fwdCover,
                    revCover = revCover,
                    metadata = list(
                      file = ifelse(is.null(bamfile),
                                     "loaded from reads",
                                     bamfile)
                      )
                   )
  }else{
    chipdata <- new("ChIPdata",
                    readsSE = emptySE,
                    readsPE = reads,
                    isPE = TRUE,
                    fwdCover = fwdCover,
                    revCover = revCover,
                    metadata = list(
                      file = ifelse(is.null(bamfile),
                                    "loaded from reads",
                                    bamfile)
                      )
                    )
    
  }
  
  chipdata

}

# 
# ##' reads class description
# ##'
# ##' Contains the reads obtained in a ChIP experiment separated by strand and then by chromosome. It has one component for each strand which are object of the data.table class with a match column to identify the regions
# ##'
# ##' @slot readsFile Character vector with the name of the reads file to be used
# ##' @slot readsF List of data.table objects containing the reads of the ChIP - Seq experiment that have + strand
# ##' @slot readsR List of data.table object containing the reads of the ChIP - Seq experiment that have - strand
# ##' @slot nReads Numeric value with the number of reads
# ##' @slot isPET Logical value indicating if the reads are paired or not
# ##'
# ##' @name reads-class
# ##' @rdname reads-class
# ##' @exportClass reads
# setClass("reads",
#   representation(
#     readsFile = "character",
#     readsF = "list",
#     readsR = "list",
#     nReads = "numeric",
#     isPET = "logical"),
#   prototype = prototype(
#     readsFile = "",
#     readsF = list(),
#     readsR = list(),
#     nReads = 0,
#     isPET = FALSE
#     ))
# 
# #setValidity("reads",
# #  function(object){
# #    return(length(object@readsF) == length(object@readsR))
# #})            
