require("GenomicRanges")
require(csaw)
#' Title
#'
#' @param file path to bam file a
#' @param file1 path to bam file b
#' @param binsize 
#' @param name sample name of file a
#' @param name1 sample name of file b
#' @param blacklist.file blacklist file in rds format
#' @param setLim plot parameter
#' @param save if True a plot is generated and saved
#' @param ext 
#'
#' @return data.frame with spearman and pearson correlation values
#' @export
#'
#' @examples
spearman_pearson_correlation <- function(file,  file1,  binsize, name, name1, blacklist.file, setLim= "auto", save=T, ext = 60){
  
  standard.chr <- paste0("chr", c(1:19, "X", "Y"))
  mm10_blacklist <-readRDS(blacklist.file)
  param <- readParam( discard=mm10_blacklist, restrict=standard.chr)
  bam.files <- c(file, file1)
  ### get nb reads in bam files
  countfile1 <- countBam(file1)
  countfile1 <- countfile1$records
  
  countfile <- countBam(file)
  countfile <- countfile$records
  ### calculate bins
  bins <- windowCounts(bam.files, bin=TRUE, width=binsize, param=param, filter=1, ext = ext)
  coverage <- assay(bins)
  coverage<- data.frame(coverage[,1]/(countfile/1000000), coverage[,2]/(countfile1/1000000))
  names(coverage) <-c("A", "B")
  ## calclate pearson and spearman correlation
  correlation.p <- cor.test(coverage[,1], coverage[,2], method = "pearson")
  corr.p.ul <- unlist(correlation.p)
  correlation.s <- cor.test(coverage[,1], coverage[,2], method = "spearman")
  corr.s.ul <- unlist(correlation.s)
  
  names.corr.p <- c(names(correlation.p)[-9], "conf.int1", "conf.int2")
  correlation.p <-data.frame(names.corr.p, corr.p.ul)
  names(correlation.p) <- c("A", "B")
  names.corr.s <- c(names(corr.s.ul)[-7])
  correlation.s <- data.frame(names.corr.s, corr.s.ul[-7])
  names(correlation.s) <- c("A", "B")
  
  ## save correlation
  correlation <- rbind(correlation.p, correlation.s)
  save.name = paste(name, "-", name1, "_correlation.norm.cvs", sep="")
  write.csv(correlation, save.name, row.names = F)
  ## save image
  ## save image
  if(save){
    if(setLim == "auto"){
      save.name = paste(name, "-", name1, "_correlation.norm.png", sep="")
      png(save.name)
      plot(coverage[,1], coverage[,2], xlab= name, ylab= name1, cex=0.7, pch = 20)
      dev.off()
    }else{
      save.name = paste(name, "-", name1, "_correlation.norm.png", sep="")
      png(save.name)
      plot(coverage[,1], coverage[,2], xlab= name, ylab= name1, xlim=c(0, setLim), 
           ylim= c(0 , setLim), cex=0.7, pch = 20)
      dev.off()
    }
  }
  
  return(correlation)
}
