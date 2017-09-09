##' @rdname isPET-methods
##' @aliases isPET
##' @docType method
##' @exportMethod isPET
setMethod("isPET",
          signature = signature(object = "ChIPdata"),
          definition = function(object)object@isPE
)

##' @rdname reads-methods
##' @aliases reads
##' @docType method
##' @exportMethod reads
setMethod("reads",
          signature = signature(object = "ChIPdata"),
          definition = function(object){
            
            if(!isPET(object)){
              object@readsSE
            }else{
              object@readsPE
            }
            
          })