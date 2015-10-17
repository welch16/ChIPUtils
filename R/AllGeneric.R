
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
##' 
##' file <- system.file("extdata","example",
##'   "encode_K562_Ctcf_first3chr_Rep1.sort.bam",package = "ChIPUtils")
##' rr <- create_reads(file)
##' readsFile(rr)
setGeneric("readsFile",
  function(object)
  standardGeneric("readsFile")
)       

##################################################################################

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
##' file <- system.file("extdata","example",
##'   "encode_K562_Ctcf_first3chr_Rep1.sort.bam",package = "ChIPUtils")
##' rr <- create_reads(file)
##' readsF(rr)
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
##' file <- system.file("extdata","example",
##'   "encode_K562_Ctcf_first3chr_Rep1.sort.bam",package = "ChIPUtils")
##' rr <- create_reads(file)
##' readsR(rr)
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

##################################################################################

##' summary methods
##'
##' summary returns a data.table object with the number of reads per strand and per chromosome
##'
##' @param object A \code{reads} object
##' @param ... Any other parameters of the summary method
##'
##' @return A data.table with number of reads by chromosome and by strand
##'
##' @export
##' @docType methods
##' @rdname summary-methods
##' @examples
##' file <- system.file("extdata","example",
##'   "encode_K562_Ctcf_first3chr_Rep1.sort.bam",package = "ChIPUtils")
##' rr <- create_reads(file)
##' summary(rr)
setGeneric("summary",
  function(object,...)
  standardGeneic("summary")
)

##################################################################################

##' nreads methods
##' 
##' nreads returns the number of reads in the object
##' 
##' @param object A \code{reads} object
##' 
##' @return The number of reads in the \code{reads} object
##' @export
##' @docType methods
##' @rdname nreads-methods
##' @examples
##' file <- system.file("extdata","example",
##'   "encode_K562_Ctcf_first3chr_Rep1.sort.bam",package = "ChIPUtils")
##' rr <- create_reads(file)
##' nreads(rr)
setGeneric("nreads",
  function(object)
  standardGeneric("nreads")
)

##################################################################################

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
##'
##' @return A \code{data.table} with two columns, shift and its respect strand cross
##' correlation
##'
##' @export
##' @docType methods
##' @rdname strand_cross_corr-methods
##' @examples
##' file <- system.file("extdata","example",
##'   "encode_K562_Ctcf_first3chr_Rep1.sort.bam",package = "ChIPUtils")
##' rr <- create_reads(file)
##' shift <- seq(10,100,by = 10)
##' chrom.sizes <- system.file("extdata","chrom.sizes","hg19.chrom.sizes",package = "ChIPUtils",mustWork = TRUE)
##' message(chrom.sizes)
##' chrom.sizes <- data.table(read.table(chrom.sizes))
##' strand_cross_corr(rr,shift, chrom.sizes)
setGeneric("strand_cross_corr",
  function(object,shift,chrom.sizes)
  standardGeneric("strand_cross_corr")
)


