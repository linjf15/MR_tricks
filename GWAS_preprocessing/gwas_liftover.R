gwas_liftover_v1 <- function(gwas,
                             col_pos_prelo="Position",
                             col_pos_poslo="Pos38",
                             col_chr="Chromosome",
                             chr_format=T,
                             hg_prelo="hg19",
                             hg_postlo="hg38")
{
  # set clock for timing
  time_start <- Sys.time()
  
  # liftOver package is required; current version is 1.18.0
  library(liftOver)
  
  # import chain file for coordination transformation
  lo_path <- system.file(package="liftOver", "extdata", paste0(hg_prelo,"To", hg_postlo ,".over.chain"))
  lo_ch <- import.chain(lo_path)
  
  # define chromosome and position columns
  pos_prelo <- gwas[[col_pos_prelo]]
  chr_prelo <- gwas[[col_chr]]
  chr_prelo <- ifelse(chr_format, paste0("chr", chr_prelo), chr_prelo)
  
  # generate GRanges object and transform
  lo_gwas_o <- GRanges(chr_prelo, IRanges(start=pos_prelo, width=1))
  seqlevelsStyle(lo_gwas_o) <- "UCSC"  # necessary
  lo_gwas_n <- liftOver(lo_gwas_o, lo_ch)
  
  # convert lo_gwas_n file as data.table object for following merging
  lo_gwas_n <- data.table(data.frame(lo_gwas_n))[,c(1,4)]
  names(lo_gwas_n) <- c("rank_id", col_pos_poslo) # define rank id
  gwas$rank_id <- seq(1:nrow(gwas)) # define rank id
  gwas <- base::merge(gwas, lo_gwas_n, by="rank_id", sort=F, all.x=T)
  
  # set clock for timing
  time_end <- Sys.time()
  cat(paste0("Time used is ", time_end-time_start," seconds\n"))
  
  # remove rank id and return gwas with liftover position
  return(gwas[,-1])
}
