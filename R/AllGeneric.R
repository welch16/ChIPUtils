
##' isPET methods
##'
##' isPET returns a logical value indicating if the object isPET
##'
##' @param object A \code{ChIPdata} object
##'
##' @return A logical value indicating if the reads are PET (TRUE) or SET (FALSE)
##'
##' @export
##' @docType methods
##' @seealso \code{\link{ChIPdata-class}}
##' @rdname isPET-methods
setGeneric("isPET",
           function(object)
             standardGeneric("isPET")
)           

##' reads methods
##' 
##' reads returns either a \code{GAlignments} or \code{GAlignmentPairs} object
##' depending on if the reads are SE or PE
##' 
##' @param object A \code{ChIPdata} object
##' 
##' @return The aligned reads of the experiment as \code{GAlignments} or 
##' \code{GAlignmentPairs} depending on if the reads are SE or PE
##' 
##' @export
##' @docType methods
##' @seealso \code{\link{ChIPdata-class}}
##' 
##' @rdname reads-methods
setGeneric("reads",
           function(object)
             standardGeneric("reads"))

##' SCC methods
##'
##' Returns a \code{tibble} with 2 columns: shift and its corresponding
##' strand cross correlation as defined in Kharchenko et al. 2008
##'
##' @param object A \code{ChIPdata} object
##' @param shift An non-negative integer vector
##' @param verbose a logical value indicating if messages should appear when calculating the SCC
##'
##' @return A \code{tibble} with two columns, shift and its respect strand cross
##' correlation
##'
##' @export
##' @docType methods
##' @rdname SCC-methods
setGeneric("SCC",
  function(object,shift,verbose)
  standardGeneric("SCC")
)

##' nreads methods
##' 
##' Returns the number of reads in the experiment
##' 
##' @param object a \code{ChIPdata} object
##' @return the number of aligned reads in the experiment
##' 
##' @export
##' @docType methods
##' @rdname nreads-methods
setGeneric("nreads",
           function(object)
             standardGeneric("nreads"))


