
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
  sizes %>% names()
}

##' @export
galignments2tibble <- function(galignment)
{
  out = as(galignment,"data.frame")
  as_tibble(out)
}

##' @export
size2GRanges <- function(tib)
{
  with(tib,
       GRanges(seqnames = chr,
               IRanges(start = 1, width = size)))
}

