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

##' @rdname nreads-methods
##' @aliases nreads
##' @docType methods
##' @exportMethod nreads
setMethod("nreads",
          signature = signature(object = "ChIPdata"),
          definition = function(object){
            reads = reads(object)
            n = length(reads)
            if(isPET(object)){
              n * 2
            }else{
              n
            }
          })