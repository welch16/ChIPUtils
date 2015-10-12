
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
##' local_strand_cross_corr(reads,region)
##' 
##' ## if the regions is empty, cross.corr = 0:
##' region <- GRanges(seqnames = "chr1",ranges = IRanges(start = 1,end = 1e4))
##' local_strand_cross_corr(reads,region,shift = 1:5)
local_strand_cross_corr <- function(reads, region,shift = 1:300 )
{
  stopifnot(length(region) == 1)
  stopifnot(class(region) == "GRanges")
  chr <- as.character(seqnames(region))
  stopifnot( chr %in% names(readsF(reads)) & chr %in% names(readsR(reads))) 

  rF <- dt2gr(readsF(reads)[[chr]])
  rR <- dt2gr(readsR(reads)[[chr]])

  ovF <- findOverlaps(region, rF)
  ovR <- findOverlaps(region, rR)
  
  if(  length(ovF)== 0 | length(ovR) == 0){
    cross.corr <- 0
  }else{
    rF <- ranges(rF[subjectHits(ovF)])
    rR <- ranges(rR[subjectHits(ovR)])

    end(rF) <- start(rF)
    start(rR) <- end(rR)

    cF <- coverage(rF,width = width(region))
    cR <- coverage(rR,width = width(region))
  
    cross.corr <- shiftApply(shift,cF,cR,cor,verbose =FALSE)
  }
  return(data.table(shift , cross.corr))
}

###############################################################################333

