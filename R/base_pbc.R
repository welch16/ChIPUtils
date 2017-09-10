
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
##' @examples
##' file <- system.file("extdata","example","encode_K562_Ctcf_first3chr_Rep1.sort.bam",package = "ChIPUtils")
##' chipdata = ChIPdata(file,isPE =  FALSE)
##'  PBC(chipdata)
PBC <- function(chipdata)
{
  stopifnot(class(chipdata) == "ChIPdata")
  alignedReads <- reads(chipdata)
  if(!isPET(chipdata)){ ## SET
    
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
        N = n()
    ) %>% ungroup()
    
    tags <- bind_rows(tagsF,tagsR)
    
  }else{ ## PET
    alignedReads <- reads(chipdata)
    alignedReads <- galignments2tibble(granges(alignedReads))
    
    tags <- alignedReads %>% 
      group_by(seqnames,start,end,strand) %>% 
      dplyr::summarise(
        N = n()
      ) %>% 
      ungroup()
    
  }
  N <- sum(tags$N)
  N1 <- tags %>% 
    filter(N == 1) %>% 
    nrow()
  N1 / N
}
