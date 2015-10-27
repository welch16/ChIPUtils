
setMethod("show",
  signature = signature(object = "reads"),
  definition = function(object){
      cat("Reads in", object@readsFile,": " ,object@nReads,"\n")
      cat("The reads are", ifelse(isPET(object), "PET", "SET"))
  }
)

###############################################################################333

##' @rdname isPET-methods
##' @aliases isPET
##' @docType method
##' @exportMethod isPET
setMethod("isPET",
  signature = signature(object = "reads"),
  definition = function(object)object@isPET
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
      stopifnot(class(value) == "list")
      object@readsF <- value
  }
)                 

###############################################################################333

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

###############################################################################333

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

###############################################################################333

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

###############################################################################333

##' @rdname strand_cross_corr-methods
##' @aliases strand_cross_corr
##' @docType methods
##' @exportMethod strand_cross_corr
setMethod("strand_cross_corr",
  signature = signature(object = "reads",shift = "numeric",
      chrom.sizes = "data.table"),
  definition = function(object,shift,chrom.sizes){

    chr <- names(readsF(object))
    dt_chr <- chrom.sizes[,1,with = FALSE]
    
    stopifnot(any(!chr %in% dt_chr))
    
    setnames(chrom.sizes,names(chrom.sizes),c("chr","size"))
    setkey(chrom.sizes,"chr")
    
    sizes <- chrom.sizes[chr]
    
    regions <- GRanges(seqnames = chr,ranges = IRanges(start = 1,
      end = sizes[,(size)]),strand = "*")
    regions <- split(regions,chr)
    scc <- lapply(regions,function(x) local_strand_cross_corr(object,x,shift) )
    sizes[,w := size / sum(size)]
    weights <- sizes[,(w)]
    scc <- mapply(function(sc,w){
      sc[,cross.corr := w * cross.corr]
      return(sc)
    },scc,weights,SIMPLIFY = FALSE)
    
    scc <- do.call(rbind,scc)
    nms <- names(scc)
    scc <- scc[,sum(cross.corr), by = shift]
    
    setnames(scc,names(scc),nms)
    return(scc)
  }
)          

###############################################################################333

