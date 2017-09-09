
galignments2tibble <- function(galignment)
{
  out = as(galignment,"data.frame")
  as_tibble(out)
}

##' Calculates the PCR bottleneck coefficient
##'
##' Calculates the PCR bottleneck coefficient defined as:
##' PBC = [# of positions with exactly 1 read mapped] / [# of positions with 1 or more reads mapped ]
##'
##' @param chipdata chipdata object
##'
##' @export
##'
##' @return Value between 0 and 1 with the PBC
##' @rdname PBC
##' @name PBC
##'
## @examples
## file <- system.file("extdata","example","encode_K562_Ctcf_first3chr_Rep1.sort.bam",package = "ChIPUtils")
## chipdata = ChIPdata(file,isPE =  FALSE)
## PBC(chipdata)
PBC <- function(chipdata)
{
  browser()
  stopifnot(class(chipdata) == "ChIPdata")
  if(!isPET(chipdata)){ ## SET
    alignedReads <- reads(chipdata)
    alignedReads <- galignments2tibble(alignedReads)
    
    tagsF <- alignedReads %>% 
      dplyr::filter(strand == "+") %>% 
      dplyr::rename(pos = start) %>% 
      group_by(seqnames,
               pos) %>% 
      dplyr::summarise(
        N = n(),
        N1 = length(unique(pos)) 
      ) %>% ungroup()
    
    tagsR <- alignedReads %>% 
      dplyr::filter(strand == "-") %>% 
      dplyr::rename(pos = end) %>% 
      group_by(seqnames,
               pos) %>% 
      dplyr::summarise(
        N = n(),
        N1 = length(unique(pos))
      ) %>% ungroup()

  }else{ ## PET
    rF <- readsF(reads)
    rR <- readsR(reads)
    rF <- lapply(rF,function(x)data.table(read = x[,paste0(start,"-",end)]))
    rR <- lapply(rR,function(x)data.table(read = x[,paste0(start,"-",end)]))
    rF <- lapply(rF,function(x)x[,V2 := 1])
    rR <- lapply(rR,function(x)x[,V2 := 1])  
    tagsF <- lapply(rF,function(x){
      out <- x[,length(V2), by = read]
      setnames(out,names(out),c("position","tags"))
      return(out)
    })
    tagsR <- lapply(rR,function(x){
      out <- x[,length(V2),by = read]
      setnames(out,names(out),c("position","tags"))
      return(out)
    })
  }
  n1 <- sum(sapply(tagsF,function(x)nrow(x[tags == 1]))) +
    sum(sapply(tagsR,function(x)nrow(x[tags == 1])))
  n2 <- sum(sapply(tagsF,nrow)) + sum(sapply(tagsR,nrow))
  return(n1 / n2)                 
}
