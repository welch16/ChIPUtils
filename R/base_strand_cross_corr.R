
##' Calculate the strand cross correlation for a specific genome region
##'
##' @description Calculates the strand cross correlation as defined by
##' Kharchenko et. al. 2008 conditional to a specific region in the genome
##'
##' @param reads A reads object
##'
##' @param region A GRanges object with length 1
##'
##' @param shift An integer vector, the default value is \code{1:300}
##'
##' @param perm A boolean value, indicating if it is gonna permute 
##' the fragments strands, the default value is \code{FALSE}
##'
##' @export
##'
##' @return A data.table with two columns \code{shift} and \code{cross.corr}
##'
##' @rdname local_strand_cross_corr
##' @name local_strand_cross_corr
##'
##' @examples
##'
##' file <- system.file("extdata","example",
##'   "encode_K562_Ctcf_first3chr_Rep1.sort.bam",package = "ChIPUtils")
##' reads <- create_reads(file)
##' region <- GRanges(seqnames = "chr1", ranges = IRanges(start = 1,end = 249250621))
##' local_strand_cross_corr(reads,region,shift = seq(10,300,by = 10))
##' 
##' ## if the regions is empty, cross.corr = 0:
##' region <- GRanges(seqnames = "chr1",ranges = IRanges(start = 1,end = 1e4))
##' local_strand_cross_corr(reads,region,shift = 1:5)
##' local_strand_cross_corr(reads,region,shift = 1:5,perm = TRUE)
local_strand_cross_corr <- function(reads, region,shift = 1:300 ,perm = FALSE)
{
  
  stopifnot(length(region) == 1)
  stopifnot(class(region) == "GRanges")
  chr <- as.character(seqnames(region))
  stopifnot( chr %in% names(readsF(reads)) & chr %in% names(readsR(reads))) 

  rF <- copy(dt2gr(readsF(reads)[[chr]]))
  rR <- copy(dt2gr(readsR(reads)[[chr]]))

  ovF <- findOverlaps(region, rF)
  ovR <- findOverlaps(region, rR)

  rF <- rF[subjectHits(ovF)]
  rR <- rR[subjectHits(ovR)]

  if(perm){
    allr <- c(rF,rR)
    strand(allr) <- sample(as.character(strand(allr)))
    allr <- split(allr,as.character(strand(allr)))
    rF <- allr[["+"]]
    rR <- allr[["-"]]
    rm(allr)
  }
  
  if(  length(ovF)== 0 | length(ovR) == 0){
    cross.corr <- 0
    shift1 <- NULL
  }else{

    end(rF) <- start(rF)
    start(rR) <- end(rR)

    rangeF <- IRanges(min(start(rF)),max(start(rF)))
    rangeR <- IRanges(min(start(rR)),max(start(rR)))
    reg <- ranges(region)

    cF <- coverage(ranges(rF))[rangeF]
    cR <- coverage(ranges(rR))[rangeR]

    ## want to make cF and cR to have the same length
    ## fix the beginning
    if(start(rangeF) != start(rangeR)){
      if(start(rangeF) < start(rangeR)){
        ext <- start(rangeR) - start(rangeF) 
        cR <- c(Rle(rep(0,ext)),cR)
      }else{
        ext <- start(rangeF) - start(rangeR) 
        cF <- c(Rle(rep(0,ext)),cF)
      }
    }
    
    ## fix the end
    if(end(rangeF) != end(rangeR)){
      if(end(rangeF) < end(rangeR)){
        ext <- end(rangeR) - end(rangeF) 
        cF <- c(cF,Rle(rep(0,ext)))
      }else{
        ext <- end(rangeF) - end(rangeR) 
        cR <- c(cR,Rle(rep(0,ext)))
      }
    }

    maxShift <- max(shift)
   
    ## fix shift
    if( maxShift >= length(cF)){
      shift1 <- shift[shift < length(cF)]    
    }else{
      shift1 <- shift
    }
    if(length(shift1) > 0){
      cc <- shiftApply(shift1,cF,cR,cor,verbose =FALSE)
    }
  }
  dt <- data.table(shift, cross.corr = 0 )
  if(length(shift1) > 0){
    dt[shift %in% shift1, cross.corr := cc]
  }
  dt[is.nan(cross.corr),cross.corr := 0]
  return(copy(dt))
}

##################################################################################

##' @rdname strand_cross_corr-methods
##' @aliases strand_cross_corr
##' @docType methods
##' @exportMethod strand_cross_corr
setMethod("strand_cross_corr",
          signature = signature(object = "reads",shift = "numeric",
                                chrom.sizes = "data.table",parallel = "logical"),
          definition = function(object,shift,chrom.sizes,parallel = FALSE){
            
            chr <- names(readsF(object))
            dt_chr <- chrom.sizes[,1,with = FALSE]
            
            if(length(chr) == 1){
              stopifnot(chr %in% as.character(dt_chr[[1]]))
            }else{
              stopifnot(any(!chr %in% dt_chr))
            }
            
            setnames(chrom.sizes,names(chrom.sizes),c("chr","size"))
            setkey(chrom.sizes,"chr")
            sizes <- chrom.sizes[chr,nomatch = 0]
            regions <- GRanges(seqnames = sizes[,(chr)],
                               ranges = IRanges(start = 1,end = sizes[,(size)]),
                               strand = "*")
            regions <- split(regions, as.character(seqnames(regions)))
            if(parallel ){
              mc <- detectCores()
              scc <- mclapply(regions,function(x) local_strand_cross_corr(object,x,shift) ,mc.cores = mc)
            }else{
              scc <- lapply(regions,function(x) local_strand_cross_corr(object,x,shift) )
            }
            weights <- copy(summary(object))
            weights[ , w := total / sum(as.numeric(total))]
            cc <- as.character(sizes[,(chr)])
            setkey(weights,chr)
            weights <- weights[chr %in%  cc ] 
            nms <- names(scc[[1]])
            scc <- do.call(rbind,scc)
            scc <- scc[,weighted.mean(cross.corr,w = weights[,(w)]),by = shift]
            setnames(scc,names(scc),nms)
            return(scc)
          }
)          

##################################################################################
  
