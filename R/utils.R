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
##' gr2 <- GRanges(seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
##'         ranges = IRanges(1:10, width = 10:1, names = head(letters,10)),
##'         strand = Rle(strand(c("-", "+", "*", "+", "-")),c(1, 2, 2, 3, 2)),
##'         score = 1:10, GC = seq(1, 0, length=10), seqinfo=seqinfo)
##' gr2dt(gr2)

gr2dt <- function(gr)
{
  out <- data.table(seqnames = as.character(seqnames(gr)),
                    start = start(gr),end = end(gr),
                    strand = as.character(strand(gr)))
  extra <- mcols(gr)
  extra <- data.table(as.data.frame(extra))
  out <- cbind(out,extra)
  return(out)
}



