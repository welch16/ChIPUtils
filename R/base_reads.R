
##' Creates a reads object
##'
##' Constructor of the reads class
##'
##' @param readsFile Name of the bam file with the aligned reads
##'
##' @export
##'
##' @return A reads object
##'
##' @rdname create_reads
##' @name create_reads
##'
##' @examples
##' file = system.file("extdata","example","encod_K562_Ctcf_first3chr_Rep1.sort.bam")
##' create_reads(file)
##' 
create_reads <- function( readsFile )
{
  bai <- paste0(readsFile,".bai")
  if(!file.exists(bai)){
    warning("Creating index file ",bai)
    message("Creating index file ",bai)
    indexBam(readsFile)
  }     
  param <- NULL
  greads <- readGAlignments(readsFile,param = param)

  greads <- as(greads,"GRanges")

  gr <- gr2dt(greads)

  seykey(gr,strand)

  fwd <- gr["+"]
  bwd <- gr["-"]

  fwd <- split(fwd,fwd[,(seqnames)])
  bwd <- split(bwd,bwd[,(seqnames)])

  out <- new("reads",readsFile = readsFile,readsF = fwd,readsR = bwd,
             nReads = nrow(gr))
  return(out)
}
