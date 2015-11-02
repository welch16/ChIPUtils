
##' Creates a reads object
##'
##' Constructor of the reads class
##'
##' @param reads_file Name of the bam file with the aligned reads
##' 
##' @param is_PET Logical indicator to asses if the reads in file are PET or SET. The default value is FALSE indicating SET reads
##'
##' @export
##'
##' @return A reads object
##'
##' @rdname create_reads
##' @name create_reads
##'
##' @examples
##' file <- system.file("extdata","example","encode_K562_Ctcf_first3chr_Rep1.sort.bam",package = "ChIPUtils",mustWork = TRUE)
##' rr <- create_reads(file)
##' summary(rr)
##' file2 <- system.file("extdata", "ex1.bam", package="Rsamtools", mustWork=TRUE)
##' create_reads(file2,is_PET = TRUE)
create_reads <- function( reads_file, is_PET = FALSE )
{
  bai <- paste0(reads_file,".bai")
  if(!file.exists(bai)){
    warning("Creating index file ",bai)
    message("Creating index file ",bai)
    indexBam(reads_file)
  }     
  param <- NULL
  if(is_PET){
    gpairs <- readGAlignmentPairs(reads_file , param = param)
    lreads <- as(left(gpairs),"GRanges")
    rreads <- as(right(gpairs),"GRanges")
    sqnms <- as.character(seqnames(lreads))
    strs <- as.character(strand(lreads))
    starts <- start(lreads)
    ends <- end(rreads)
    idx <- which(ends - starts + 1 > 0)
    stopifnot(length(idx) > 0)
    greads <- GRanges(seqnames = sqnms[idx],
      ranges = IRanges(start = starts[idx],end = ends[idx]),
      strand = strs[idx])
  }else{
    greads <- readGAlignments(reads_file,param = param)
    greads <- as(greads,"GRanges")
  }
  
  gr <- gr2dt(greads)
  setkey(gr,strand)

  fwd <- gr["+",nomatch = 0]
  bwd <- gr["-",nomatch = 0]
  fwd <- split(fwd,fwd[,(seqnames)])
  bwd <- split(bwd,bwd[,(seqnames)])

  out <- new("reads",readsFile = reads_file,readsF = fwd,readsR = bwd,
             nReads = nrow(gr),isPET = is_PET)
  return(out)
}
