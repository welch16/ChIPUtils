##' @importFrom methods setClass setGeneric setMethod setRefClass
NULL

##' reads class description
##'
##' Contains the reads obtained in a ChIP experiment separated by strand and then by chromosome. It has one component for each strand which are object of the data.table class with a match column to identify the regions
##'
##' @slot readsFile Character vector with the name of the reads file to be used
##' @slot readsF List of data.table objects containing the reads of the ChIP - Seq experiment that have + strand.
##' @slot readsR List of data.table object containning the reads of the ChIP - Seq experiment that have - strand
##' @slot nReads Numeric value with the number of reads
##' @slot isPET Logical variable indicating if the reads the number of sequenced ends in the fragment (TRUE = Paired, FALSE = Single)
##'
##' @name reads-class
##' @rdname reads-class
##' @exportClass reads
setClass("reads",
  representation(
    readsFile = "character",
    readsF = "list",
    readsR = "list",
    nReads = "numeric",
    isPET = "logical"),
  prototype = prototype(
    readsFile = "",
    readsF = list(),
    readsR = list(),
    nReads = 0,
    isPET = FALSE
    ))

setValidity("reads",
  function(object){
    return(length(object@readsF) == length(object@readsR))
})            
