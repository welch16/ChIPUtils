##' @import data.table
##' @importFrom GenomicRanges GRanges
##' @importFrom IRanges IRanges
NULL

##' Convert a GRanges object to data.table
##'
##' @description Converts a GRanges object to data.table, the columns
##' are gonna have the same names as the GRanges fields
##'
##' @param gr A GRanges object
##'
##' @export
##'
##' @return A data.table object
##'
##' @rdname gr2dt
##' @name gr2dt
##'
##' @seealso \code{\link{dt2gr}}
##' @examples
##'
##' gr <- GRanges(seqnames = "chr1",
##'         ranges = IRanges(start = 3:10,width = 30),
##'         strand = "*")
##' gr2dt(gr)
##'

gr2dt <- function(gr)
{
  out <- data.table(seqnames = as.character(seqnames(gr)),
                    start = start(gr),end = end(gr),
                    strand = as.character(strand(gr)))
  extra <- mcols(gr)
  extra <- data.table(as.data.frame(extra))
  if(nrow(extra)>0){
    out <- cbind(out,extra)
  }
  return(out)
}



