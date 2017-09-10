
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
##' @seealso \code{\link{ChIPData-class}}
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
##' @seealso \code{\link{ChIPData-class}}
##' 
##' @rdname reads-methods
setGeneric("reads",
           function(object)
             standardGeneric("reads"))

##' strand_cross_corr methods
##'
##' Returns a \code{data.table} with 2 columns: shift and its corresponding
##' strand cross correlation as defined in Kharchenko et al. 2008
##'
##' @param object A \code{reads} object
##' @param shift An non-negative integer vector
##' @param chrom.sizes An element representing the sizes of the chromosomes for the
##' genome of interest. It could be a \code{data.table} or a\code{character} pointing
##' to any of the genomes that are available in the package. The \code{data.table} must have two columns
##' the chromosome and its size.
##' @param parallel A logical value indicating if the process is gonna be run in parallel
##'
##' @return A \code{data.table} with two columns, shift and its respect strand cross
##' correlation
##'
##' @export
##' @docType methods
##' @rdname strand_cross_corr-methods
##' @examples
##' file <- system.file("extdata","example","encode_K562_Ctcf_first3chr_Rep1.sort.bam",package = "ChIPUtils")
##' rr <- create_reads(file)
##' shift <- seq(10,100,by = 10)
##' chrom.sizes <- system.file("extdata","chrom.sizes","hg19.chrom.sizes",package = "ChIPUtils",mustWork = TRUE)
##' chrom.sizes <- data.table(read.table(chrom.sizes))
##' strand_cross_corr(rr,shift, chrom.sizes,FALSE)
setGeneric("strand_cross_corr",
  function(object,shift,chrom.sizes,parallel)
  standardGeneric("strand_cross_corr")
)


