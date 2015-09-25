##' @importFrom methods setClass setGeneric setMethod setRefClass
NULL

##' reads class description
##'
##' Contains the reads obtained in a ChIP - seq experiment separated by strand and then by chromosome. It has one component for each strand which are object of the data.table class with a match column to identify the regions
##'
##' @slot readsF List of data.table objects containing the reads of the ChIP - Seq experiment that have + strand.
##' @slot readsR List of data.table object containning the reads of the ChIP - Seq experiment that have - strand
##' @seealso \code{\link{loadReads}}
##'
##' @name reads-class
##' @rdname reads-class
##' @exportClass reads
setClass("reads",
  representation(readsF = "list",readsR = "list"),
  prototype = prototype(readsF = list(),readsR = list()))

setValidity("reads",
  function(object){
    return(length(object@readsF) == length(object@readsR))
})            
