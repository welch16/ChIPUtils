
##' Calculates the PCR bottleneck coefficient
##'
##' Calculates the PCR bottleneck coefficient defined as:
##' PBC = [# of positions with exactly 1 read mapped] / [# of positions with 1 or more reads mapped ]
##'
##' @param reads Reads object
##'
##' @export
##'
##' @return Value between 0 and 1 with the PBC
##' @rdname PBC
##' @name PBC
##'
##' @examples
##' file <- system.file("extdata","example",
##'   "encode_K562_Ctcf_first3chr_Rep1.sort.bam",package = "ChIPUtils")
##' rr <- create_reads(file)
##' PBC(rr)

PBC <- function(reads)
{
  stopifnot(class(reads) == "reads")
  tagsF <- lapply(readsF(reads),
    function(x){
      out <- x[,length(end),by = start]
      setnames(out,names(out),c("position","tags"))
      return(out)})
  tagsR <- lapply(readsR(reads),
    function(x){
      out <- x[,length(start),by = end]
      setnames(out,names(out),c("position","tags"))
      return(out)})
  n1 <- sum(sapply(tagsF,function(x)nrow(x[tags == 1]))) +
      sum(sapply(tagsR,function(x)nrow(x[tags == 1])))
  n2 <- sum(sapply(tagsF,nrow)) + sum(sapply(tagsR,nrow))
  return(n1 / n2)                 
}
