
##' Returns the chromosome sizes available in the package
##'
##' @description Returns the chromosome sizes stored in the \code{inst/extdata/chrom.size} directory of the ChIPUtils package
##'
##' @export
##'
##' @return A character vector
##'
##' @rdname possibleChrSizes
##' @name possibleChrSizes
##'
##' @examples
##'
##' possibleChrSizes()
possibleChrSizes <- function()
{
  data(chromSizes)
  chromSizes %>% names()
}

galignments2tibble <- function(galignment)
{
  out = as(galignment,"data.frame")
  as_tibble(out)
}

size2GRanges <- function(tib)
{
  with(tib,
       GRanges(seqnames = chr,
               IRanges(start = 1, width = size)))
}

genomeFromChIPdata <- function(chipdata)
{
  reads <- granges(reads(chipdata))
  chr <- seqlevelsInUse(reads)
  starts <- rep_len(1,length(chr))
  GRanges(
    seqnames = chr,
    IRanges(starts,seqlengths(reads))
  )
  
}
