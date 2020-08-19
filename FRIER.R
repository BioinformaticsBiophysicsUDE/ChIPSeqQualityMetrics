require("csaw")
require("GenomicAlignments")
require("foreach")

#' Title
#'
#' @param bam path to bam file
#' @param blacklist file with blacklisted regions in rds-format
#' @param minq minimum mapping quality
#' @param window.width 
#' @param bin.width 
#' @param filter 
#' @param fraglen 
#'
#' @return FRIER data frame
#' @export
#'
#' @examples
calculateEnrichedRanges <- function(bam, 
                                    blacklist, minq=25,
                                    window.width=150, bin.width=3000, filter=1, fraglen=60){

  # initialize parameter
  standard.chr <- paste0("chr", c(1:19, "X", "Y"))
  blacklist <- readRDS(blacklist)
  param <- readParam(minq=minq, discard=blacklist, restrict=standard.chr)
  ## if no fraglen is provided, fraglength is calculated
  if(length(fraglen) == 0){
    x <- correlateReads(bam, param=reform(param, dedup=TRUE))
    fraglen <- which.max(x) - 1
    cat("calculated fragment length: ", fraglen, "\n")
  }
  ## get enriched regions
  win.data <- windowCounts(bam, param=param, width=window.width, ext=fraglen, filter=filter)
  #Filtering windows by abundance
  # Count the number of extended reads overlapping a sliding window at spaced positions across the genome.
  bins <- windowCounts(bam, bin=TRUE, width=bin.width, param=param) 
  # for higher R version
  # filter.stat <- filterWindowsGlobal(win.data, bins)
  filter.stat <- filterWindows(win.data, bins, type="global")
  bam2 <- readGAlignments(bam)
  min.fc <- c(seq(1.1 ,16.1, by = 0.5))
  
  FRIER=foreach(i =  min.fc, .combine = "c") %do% {
    keep <- filter.stat$filter > log2(i)
    summay.keep = summary(keep)
    #print(summay.keep)
    #The actual filtering itself is done by simply subsetting the RangedSummarizedExperiment object.
    filtered.data <- win.data[keep,]
    retained.ranges <- filtered.data@rowRanges
    retained.ranges.reduced <- reduce(retained.ranges)
    #summarizeOverlaps -  simplifies counting reads in genomic ranges 
    summary.overlap <- data.frame(intNotEmpty = 
                                    assay(summarizeOverlaps(GRangesList(retained.ranges.reduced), bam2, 
                                                            mode="IntersectionNotEmpty"), minoverlap=1, type = "any"))
    summary.overlap$reads/length(bam2)
  }
  FRIER.m <- as.data.frame(t(matrix(c(log2(min.fc),FRIER), nrow=2, byrow = T)))
  names(FRIER.m) <- c("fold_enrichment", "reads_%")
  return(FRIER.m)
}
