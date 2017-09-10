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

# 
# setMethod("show",
#   signature = signature(object = "reads"),
#   definition = function(object){
#       cat("Reads in", object@readsFile,": " ,object@nReads,"\n")
#       cat("The reads are", ifelse(isPET(object), "PET", "SET\n"))
#   }
# )
# 
# 
# ##################################################################################
# 
# ##' @rdname summary-methods
# ##' @aliases summary
# ##' @docType methods
# ##' @exportMethod summary
# setMethod("summary",
#   signature = signature(object = "reads"),
#   definition = function(object){
#     nreadsF <- sapply(readsF(object),nrow)
#     nreadsR <- sapply(readsR(object),nrow)
#     out <- data.table(chr = names(readsF(object)),readsF = nreadsF,readsR = nreadsR,
#       total = nreadsF + nreadsF)
#     out
#   }
# )
# 
