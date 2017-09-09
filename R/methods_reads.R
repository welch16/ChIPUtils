
setMethod("show",
  signature = signature(object = "reads"),
  definition = function(object){
      cat("Reads in", object@readsFile,": " ,object@nReads,"\n")
      cat("The reads are", ifelse(isPET(object), "PET", "SET\n"))
  }
)

##################################################################################

##' @rdname readsF-methods
##' @aliases readsF
##' @docType methods
##' @exportMethod readsF
setMethod("readsF",
  signature = signature(object = "reads"),
  definition = function(object)object@readsF
)

##' @rdname readsF-methods
##' @aliases readsF<-
##' @docType methods
##' @exportMethod readsF<-
setReplaceMethod("readsF",
  signature = signature(object = "reads",value = "list"),
  definition = function(object,value){
      stopifnot(class(value) == "list")
      object@readsF <- value
  }
)                 

##################################################################################

##' @rdname readsR-methods
##' @aliases readsR
##' @docType methods
##' @exportMethod readsR
setMethod("readsR",
  signature = signature(object = "reads"),
  definition = function(object)object@readsR
)

##' @rdname readsR-methods
##' @aliases readsR<-
##' @docType methods
##' @exportMethod readsR<-
setReplaceMethod("readsR",
  signature = signature(object = "reads",value = "list"),
  definition = function(object,value){
      stopifnot(class(value) == "list")
      object@readsR <- value
  }
)                 

##################################################################################

##' @rdname summary-methods
##' @aliases summary
##' @docType methods
##' @exportMethod summary
setMethod("summary",
  signature = signature(object = "reads"),
  definition = function(object){
    nreadsF <- sapply(readsF(object),nrow)
    nreadsR <- sapply(readsR(object),nrow)
    out <- data.table(chr = names(readsF(object)),readsF = nreadsF,readsR = nreadsR,
      total = nreadsF + nreadsF)
    out
  }
)

##################################################################################

##' @rdname nreads-methods
##' @aliases nreads
##' @docType methods
##' @exportMethod nreads
setMethod("nreads",
  signature = signature(object = "reads"),
  definition = function(object){
    return(object@nReads)
  }
)

##################################################################################



