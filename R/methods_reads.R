
setMethod("show",
  signature = signature(object = "reads"),
  definition = function(object){
      cat("Reads in ", object@readsFile,"\n")
      cat("It contains ", object@nReads,"aligned reads\n")
  }
)

###############################################################################333

##' @rdname readsFile-methods
##' @aliases readsFile
##' @docType method
##' @exportMethod readsFile
setMethod("readsFile",
  signature = signature(object = "reads"),
  definition = function(object)object@readsFile
)

###############################################################################333

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
      stopifnot(class(object) == "list")
      object@readsF <- value
  }
)                 

###############################################################################333

##' @rdname readsR-methods
##' @aliases readsR
##' @docType methods
##' @exportMethod readsR
setMethod("readsF",
  signature = signature(object = "reads"),
  definition = function(object)object@readsR
)

##' @rdname readsR-methods
##' @aliases readsR<-
##' @docType methods
##' @exportMethod readsF<-
setReplaceMethod("readsF",
  signature = signature(object = "reads",value = "list"),
  definition = function(object,value){
      stopifnot(class(object) == "list")
      object@readsR <- value
  }
)                 
