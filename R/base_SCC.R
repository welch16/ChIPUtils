
##' @rdname SCC-methods
##' @aliases SCC
##' @docType methods
##' @exportMethod SCC
setMethod("SCC",
          signature = signature(object = "ChIPdata",shift = "numeric",
                                verbose = "logical"),
          definition = function(object,shift,verbose = FALSE){
            
            mc = getOption("mc.cores")
            if(is.null(mc)){
              mc = 1
            }
            
            if(Sys.info()[["sysname"]] == "Windows"){
              myparam = SnowParam(workers = mc,type = "SOCK")
            }else{
              myparam = MulticoreParam(workers = mc)
            }
            
            scc_list <- bpmapply(function(fwdCover,revCover,nm,shift){
              if(verbose)message(nm,":")
              shiftApply(shift,fwdCover,revCover,cor,verbose =verbose)
            },
            object@fwdCover,
            object@revCover,
            names(object@fwdCover),
            MoreArgs = list(shift),SIMPLIFY = FALSE,
            BPPARAM = myparam)
            
            chrs <- names(object@fwdCover)
            scc <- mapply(function(x,y)tibble(shift,corr = x,chr = y),
                          scc_list,
                          chrs,SIMPLIFY = FALSE) %>% 
              bind_rows()
            
            if(isPET(object)){
              weights <- galignments2tibble(granges(reads(object)))
            }else{
              weights <- galignments2tibble(reads(object))
            }
            weights <- weights %>% 
              dplyr::group_by(seqnames) %>% 
              dplyr::summarise(
                nreads = n()
              ) %>% 
              dplyr::mutate(
                weights = nreads / sum(nreads)
              ) %>% 
              dplyr::rename(chr =seqnames) %>% 
              ungroup()
            
            w <- weights$weights

            scc <- scc %>% 
              group_by(shift) %>% 
              summarise(
                corr = weighted.mean(corr,w)
              ) %>% 
              ungroup()
            scc
          }
)          
