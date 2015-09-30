
##' readsFile methods
##'
##' readsFile returns the name of the file with the aligned reads
##'
##' @param object A \code{reads} object
##'
##' @return A character with the name of the file with the aligned reads
##'
##' @export
##' @docType methods
##' @seealso \code{\link{reads-class}}
##' @rdname readsFile-methods
##' @examples
##' \dontrun{
##' readsFile(reads)
##' }
setGeneric("readsFile",
  function(object)
  standardGeneric("readsFile")
)       

###############################################################################333

##' readsF methods
##'
##' readsF returns a list of the fowward aligned reads 
##'
##' @param object A \code{reads} object
##'
##' @return A list of data.table with same format as a GenomicRanges object with an additional match column
##'
##' @export
##' @docType methods
##' @seealso \code{\link{reads-class}}
##' @rdname readsF-methods
##' @examples
##' \dontrun{
##' readsF(segvis)
##' readsF(reads)
##' readsF(segvis) <- new_reads
##' }
setGeneric("readsF",
  function(object)
  standardGeneric("readsF")
)           

##' readsF<- assisgn a list of forward reads to a reads object
##'
##' @param value A list of data.table with same format as a GenomicRanges object with an additional match column
##'
##' @return A reads object with modified forward reads
##' @rdname readsF-methods
setGeneric("readsF<-",
  function(object,value)
  standardGeneric("readsF<-")
)           

###############################################################################333

##' readsR methods
##'
##' readsR returns a list of the backward aligned reads 
##'
##' @param object A \code{reads} object
##'
##' @return A list of data.table with same format as a GenomicRanges object with an additional match column
##'
##' @export
##' @docType methods
##' @seealso \code{\link{reads-class}}
##' @rdname readsR-methods
##' @examples
##' \dontrun{
##' readsF(segvis)
##' readsF(reads)
##' readsF(segvis) <- new_reads
##' }
setGeneric("readsR",
  function(object)
  standardGeneric("readsR")
)           

##' readsR<- assisgn a list of forward reads to a reads object
##'
##' @param value A list of data.table with same format as a GenomicRanges object with an additional match column
##'
##' @return A reads object with modified backward reads
##' @rdname readsR-methods
setGeneric("readsR<-",
  function(object,value)
  standardGeneric("readsR<-")
)           

