##' Calculates the Standarized Standard Deviation
##'
##' Calculates the SSD coefficient as defined in Planet et. al. 2011 as a uniformity measure
##'
##' @param reads Reads object
##'
##' @export
##'
##' @return A non-negative value
##' @rdname SSD
##' @name SSD
##'
##' @examples
##' file <- system.file("extdata","example","encode_K562_Ctcf_first3chr_Rep1.sort.bam",package = "ChIPUtils")
##' rr <- create_reads(file)
##' SSD(rr)
SSD <- function(reads)
{
  stopifnot(class(reads) == "reads")
  rF <- lapply(readsF(reads),dt2ir)
  rR <- lapply(readsR(reads),dt2ir)
  gr <- mapply(c,rF,rF,SIMPLIFY = FALSE)
  cover <- lapply(gr,coverage)
  w <- sapply(cover,length)
  out <- weighted.mean(sapply(cover,sd),w=w)
  out <- out / sqrt(sum(w))
  return(out)
}
