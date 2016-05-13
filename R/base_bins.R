##' @importFrom scales trans_format
##' @importFrom scales math_format
NULL

##' Partitions a genome into bins of a fixed length and counts the number of reads that overlap each bin
##'
##' @description Partitions a chromosome into fixed length bins and counts the number of extended fragments that 
##' overlap each bin
##'
##' @param bin_size A integer value used to partition the chromosome 
##'
##' @param reads A reads object
##'
##' @param chrom A GRanges object specifying the genome to bin. The maximum length used to create the bins
##' is gonna be used the integer part of (chromLen / fragLen) times fragLen for each chromosome.
##' 
##' @param frag_len An integer value used to extend the fragments. The default value is one, to count only the 5' ends 
##' that overlaps the bins
##'
##' @export
##'
##' @return A GRanges object. If a reads object is present then is going to count the number of extended fragments per 
##' bin
##'
##' @rdname create_bins
##' @name create_bins
##'
##' @examples
##'
##' file <- system.file("extdata","example",
##'   "encode_K562_Ctcf_first3chr_Rep1.sort.bam",package = "ChIPUtils")
##' reads <- create_reads(file)
##' chrom <- GRanges(seqnames = "chr1", ranges = IRanges(start = 1,end = 249250621))
##' create_bins(2000,reads = reads)
##' create_bins(2000, chrom = chrom )
##' create_bins(2000, reads = reads, frag_len = 200)
##' create_bins(2000, reads = reads,chrom = chrom, frag_len = 200)
create_bins <- function(bin_size, reads = NULL , chrom = NULL, frag_len = 1)
{
  stopifnot( !is.null(reads) | !is.null(chrom))
  stopifnot( frag_len > 0)
  if(frag_len > 1 & is.null(reads)){
    message("The reads are absent, frag_len is going to be ignored")
  }
  stopifnot( bin_size > 0)
  if(!is.null(reads)){
    chr <- names(readsF(reads))
    all_reads <- mapply(rbind,readsF(reads),readsR(reads),SIMPLIFY = FALSE)
    if(is.null(chrom)){
      starts <- sapply(all_reads,function(x) min(x[,min(start)], x[,min(end)] ))
      ends <- sapply(all_reads,function(x) max(x[,max(end)] , x[,max(start)]))
      chrom <- GRanges(chr , ranges = IRanges(start = starts , end = ends),strand = "*")
      rm(starts,ends)
    }
  }
  if(!is.null(chrom)){
    stopifnot(class(chrom) == "GRanges")
    chr <- as.character(seqnames(chrom))
  }

  chrom <- split(chrom,as.character(seqnames(chrom)))
  coord <- lapply(chrom,function(x)
    seq( trunc(start(x) / bin_size) * bin_size, trunc(end(x)/bin_size) * bin_size , by = bin_size))
  
  bins <- mapply(function(chr,starts,bin_size)
    GRanges(seqnames = chr, ranges = IRanges(start = starts , width = bin_size),strand = "*" )
  ,chr,coord, MoreArgs = list(bin_size),SIMPLIFY = FALSE)  
  
  rm(coord)  
  names(bins) <- NULL
  suppressWarnings(bins <- do.call(c,bins))
  
  if(!is.null(reads)){
    gr <- lapply(all_reads,dt2gr)
    names(gr) <- NULL
    suppressWarnings(gr <- do.call(c,gr))
    if(!isPET(reads)){
      if(frag_len > 1){
        gr <- resize(gr,frag_len)
      }
    }else{
      message("The reads are PET. Therefore the fragments shouldn't be extended")
    }
    tagCounts <- countOverlaps(bins,gr)
    mcols(bins)[["tagCounts"]] <- tagCounts
  }

  return(bins)
  
}

###############################################################################333

##' Creates a quick hexbin plot to check the relationship between two reads objects
##' 
##' @description Creates a quick hexbin plot to visualize the relationship between 
##' two reads objects
##' 
##' @param reads_x A reads object that corresponds to the x-axis
##' 
##' @param reads_y Another reads object, this one correspond to the y-axis
##' 
##' @param bin_size A integer value used to partition the chromosome 
##' 
##' @param log A logical flag that indicates if log10 scale is going to be used in the axes. The default value is FALSE
##' 
##' @param ma A logical flag that indicates if there is going to be the change of variable M := log2(x*y) and A := log2(x/y),
##' in the case that both log and ma are true, log is going to be ignored
##' 
##' @param nr_bins Integer value with the number of bins used to build the plot. The default value is 100
##'
##' @param chrom A GRanges object specifying the genome to bin. The maximum length used to create the bins
##' is gonna be used the integer part of (chromLen / fragLen) times fragLen for each chromosome.
##' 
##' @param frag_len An integer value used to extend the fragments. The default value is one, to count only the 5' ends 
##' that overlaps the bins
##'
##' @export
##' 
##' @return A ggplot object with the hexbin plot
##' 
##' @rdname hexbin_plot
##' @name hexbin_plot
##' 
##' @examples 
##' \dontrun{
##' file_x <- system.file("extdata","example",
##'   "encode_K562_Ctcf_first3chr.sort.bam",package = "ChIPUtils")
##' file_y <- system.file("extdata","example",
##'   "encode_K562_H3k27ac_first3chr.sort.bam",package = "ChIPUtils")
##'   
##' reads_x <- create_reads(file_x)
##' reads_y <- create_reads(file_y)
##'
##' hexbin_plot(reads_x,reads_y,1e3)+xlim(0,500)+ylim(0,500)
##' hexbin_plot(reads_x,reads_y,1e3,frag_len = 2000)+xlim(0,500)+ylim(0,500)
##' hexbin_plot(reads_x,reads_y,1e3,frag_len = 2000,log = TRUE)
##' hexbin_plot(reads_x,reads_y,1e3,frag_len = 2000,ma = TRUE)
##' }
hexbin_plot <- function(reads_x,reads_y,bin_size,log = FALSE,ma = FALSE,nr_bins = 100,chrom = NULL, frag_len = 1)
{
  stopifnot(class(reads_x) == "reads")
  stopifnot(class(reads_y) == "reads")
  stopifnot(bin_size > 0)
  stopifnot(frag_len > 0)
  stopifnot(nr_bins > 1)
  
  if(is.null(chrom)){

    chr_x <- intersect(
      names(readsF(reads_x)),names(readsR(reads_x)))
    chr_y <- intersect(
      names(readsF(reads_y)),names(readsR(reads_y)))
    common <- intersect(chr_x,chr_y)

    all_reads <- mapply(rbind,
      readsF(reads_x)[common],readsR(reads_x)[common],
      readsF(reads_y)[common],readsR(reads_y)[common],SIMPLIFY = FALSE)
    starts <- sapply(all_reads,
      function(x) min(x[,min(start)], x[,min(end)] ))
    ends <- sapply(all_reads,function(x) max(x[,max(end)] , x[,max(start)]))
    chr <- names(starts)
    
    chrom <- GRanges(chr ,
      ranges = IRanges(start = starts , end = ends),strand = "*")
    rm(starts,ends,all_reads,chr,chr_x,chr_y,common)
    
  }

  bins_x <- create_bins(bin_size,reads_x,chrom = chrom,  frag_len = frag_len )
  bins_y <- create_bins(bin_size,reads_y,chrom = chrom,  frag_len = frag_len )
  
  dt <- data.table(x = mcols(bins_x)$tagCounts,y = mcols(bins_y)$tagCounts)
  
  if(ma & log){
    message("Both ma and log are TRUE, log is going to be ignore")
  }
  
  if(log & !ma){
    dt[,x := 1 + x]
    dt[,y := 1 + y]
  }
  
  if(ma){
    dt <- dt[ x > 0 & y > 0]
    dt[,M := log2(x*y)]
    dt[,A := log2(x/y)]
    dt[,x := NULL]
    dt[,y := NULL]
    setnames(dt,names(dt),c("x","y"))
  }
  
  r <- viridis::viridis(100, option = "D")
  p <- ggplot(dt, aes(x,y))+stat_binhex(bins = nr_bins)+
    scale_fill_gradientn(colours = r,trans = 'log10',
      labels=trans_format('log10',math_format(10^.x)) )
  
  if(log & !ma){
    p <- p + scale_x_log10()+scale_y_log10()
  }  
  
  if(ma){
    p <- p + xlab("M") + ylab("A")
  }
  
  return(p)
}










