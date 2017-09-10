##' @importFrom scales trans_format
##' @importFrom scales math_format
NULL

##' Partitions a genome into bins of a fixed length and counts the number of reads that overlap each bin
##'
##' @description Partitions a chromosome into fixed length bins and counts the number of extended fragments that 
##' overlap each bin
##'
##' @param binSize A integer value used to partition the chromosome
##' @param chipdata A \code{ChIPdata} object
##' @param chrom A GRanges object specifying the genome to bin. The maximum length used to create the bins
##' is gonna be used the integer part of (chromLen / fragLen) times fragLen for each chromosome.
##' 
##' @param fragLen An integer value used to extend the fragments. The default value is one, to count only the 5' ends 
##' that overlaps the bins. For PE reads is going to be ignored
##'
##' @export
##' @return A GRanges object. If a reads object is present then is going to count the number of extended fragments per 
##' bin
##'
##' @rdname createBins
##' @name createBins
createBins <- function(binSize, chipdata = NULL , chrom = NULL, fragLen = 1)
{
  stopifnot( !is.null(chipdata) | !is.null(chrom))
  stopifnot( fragLen > 0)
  if(fragLen > 1 & is.null(chipdata)){
    warning("The reads are absent, frag_len is going to be ignored")
  }
  stopifnot( binSize > 0)
  
  if(!is.null(chipdata)){
    reads <- reads(chipdata)
    if(isPET(chipdata)){
      reads <- granges(reads)
    }
    chr <- seqlevelsInUse(reads)
    if(is.null(chrom)){
      chrom <- genomeFromChIPdata(chipdata)
    }
  }
  if(!is.null(chrom)){
    stopifnot(class(chrom) == "GRanges")
    chr <- seqlevelsInUse(chrom)
    seqlevels(chrom) <- chr
  }

  chrom <- S4Vectors::split(chrom,seqnames(chrom))
  coord <- lapply(chrom,function(x){
    seq( trunc(start(x) / binSize) * binSize,
         trunc(end(x)/binSize) * binSize , by = binSize)})
  
  bins <- mapply(function(chr,starts,binSize)
    GRanges(seqnames = chr, 
            IRanges(start = starts , width = binSize),
            strand = "*" ),
    chr,coord, MoreArgs = list(binSize),SIMPLIFY = FALSE)  
  
  rm(coord)  
  names(bins) <- NULL
  bins <- unlist(GRangesList(bins))

  if(!is.null(chipdata)){
    gr <- granges(reads(chipdata))
    if(!isPET(chipdata)){
      if(fragLen > 1){
        gr <- resize(gr,fragLen)
      }
    }else{
      warning("The reads are PET. Therefore the are not going to be extended")
    }
    mcols(bins)[["tagCounts"]] <- countOverlaps(bins,gr)
  }

  bins

}


##' Creates a quick hexbin plot to check the relationship between two ChIPdata objects
##' 
##' @description Creates a quick hexbin plot to visualize the relationship between 
##' two ChIPdata objects
##' 
##' @param x A \code{ChIPdata} object that corresponds to the x-axis
##' @param y Another \code{ChIPdata} object, this one correspond to the y-axis
##' @param binSize A integer value used to partition the chromosome 
##' @param log A logical flag that indicates if log10 scale is going to be used in the axes. The default value is FALSE
##' @param nrBins Integer value with the number of bins used to build the plot. The default value is 100
##' @param chrom A GRanges object specifying the genome to bin. The maximum length used to create the bins
##' is gonna be used the integer part of (chromLen / fragLen) times fragLen for each chromosome.
##' @param fragLen An integer value used to extend the fragments. The default value is one, to count only the 5' ends 
##' that overlaps the bins
##'
##' @export
##' 
##' @return A ggplot object with the hexbin plot
##' 
##' @rdname hexbinPlot
##' @name hexbinPlot
hexbinPlot <- function(x,y,binSize,log = FALSE,nrBins =40,chrom = NULL, fragLen = 1)
{
  stopifnot(class(x) == "ChIPdata",
            class(y) == "ChIPdata")
  stopifnot(binSize > 0,
            fragLen > 0,
            nrBins > 1)

  if(is.null(chrom)){
  
    chromx <- genomeFromChIPdata(x)
    chromy <- genomeFromChIPdata(y)
    
    chrom <- GenomicRanges::intersect(chromx,chromy)
  
  }

  binsx <- createBins(binSize,x,chrom = chrom,fragLen = fragLen)
  binsy <- createBins(binSize,y,chrom = chrom,fragLen = fragLen)
  
  plotData <- tibble(
    x = mcols(binsx)$tagCounts,
    y = mcols(binsy)$tagCounts
  )

  if(log){
    
    plotData <- plotData %>% 
      mutate_all(funs(1 + x))
    
  }

  pal <- viridis::viridis(100, option = "D")
  plot <- plotData %>% 
    ggplot(aes(x,y))+
    stat_binhex(bins = nrBins)+
    scale_fill_gradientn(colours = pal,trans = 'log10',
      labels=trans_format('log10',math_format(10^.x)) )
  
  if(log){
    plot <- plot + scale_x_log10()+scale_y_log10()
  }  
  
  plot
}

##' Creates a quick MA plot to check the relationship between two reads objects
##' 
##' @description Creates a quick MA plot to visualize the relationship between 
##' two reads objects
##' 
##' @param x A \code{ChIPdata} object that corresponds to the x-axis
##' @param y Another \code{ChIPdata} object, this one correspond to the y-axis
##' @param binSize A integer value used to partition the chromosome 
##' @param nrBins Integer value with the number of bins used to build the plot. The default value is 100
##' @param chrom A GRanges object specifying the genome to bin. The maximum length used to create the bins
##' is gonna be used the integer part of (chromLen / fragLen) times fragLen for each chromosome.
##' @param fragLen An integer value used to extend the fragments. The default value is one, to count only the 5' ends 
##' that overlaps the bins
##'
##' @export
##' 
##' @return A ggplot object with the hexbin plot
##' 
##' @rdname MAPlot
##' @name MAPlot
MAPlot <- function(x,y,binSize,nrBins =40,chrom = NULL, fragLen = 1)
{
  stopifnot(class(x) == "ChIPdata",
            class(y) == "ChIPdata")
  stopifnot(binSize > 0,
            fragLen > 0,
            nrBins > 1)
  
  if(is.null(chrom)){
    
    chromx <- genomeFromChIPdata(x)
    chromy <- genomeFromChIPdata(y)
    
    chrom <- GenomicRanges::intersect(chromx,chromy)
    
  }
  
  binsx <- createBins(binSize,x,chrom = chrom,fragLen = fragLen)
  binsy <- createBins(binSize,y,chrom = chrom,fragLen = fragLen)
  
  plotData <- tibble(
    x = mcols(binsx)$tagCounts,
    y = mcols(binsy)$tagCounts
  ) %>% 
    mutate(
      M = 0.5 * (log2(x) + log2(y)),
      A = log2(x) - log2(y)
    )
  
  

  pal <- viridis::viridis(100, option = "D")
  plot <- plotData %>% 
    ggplot(aes(M,A))+
    stat_binhex(bins = nrBins)+
    scale_fill_gradientn(colours = pal,trans = 'log10',
                         labels=trans_format('log10',math_format(10^.x)) )
  
  plot
}









