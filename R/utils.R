
##' Returns the chromosome sizes available in the package
##'
##' @description Returns the chromosome sizes stored in the \code{inst/extdata/chrom.size} directory of the ChIPUtils package
##'
##' @export
##'
##' @return A character vector
##'
##' @rdname possible_chrom_sizes
##' @name possible_chrom_sizes
##'
##' @examples
##'
##' possible_chrom_sizes()

possible_chrom_sizes <- function()
{
  dr <- system.file("inst","extdata","chrom.sizes",package = "ChIPUtils")
  sizes <- list.files(dr)
  chrom <- sapply(strsplit(sizes,".",fixed = TRUE),function(x)x[1])
  return(chrom)

}

galignments2tibble <- function(galignment)
{
  out = as(galignment,"data.frame")
  as_tibble(out)
}
