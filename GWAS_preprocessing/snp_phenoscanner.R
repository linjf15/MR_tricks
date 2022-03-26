# built for quickly queries the PhenoScanner database
# the input dat can be either a data.frame containing SNP column
# or just SNP list
  
snp_phenoscanner <- function(dat = exposure_MAP)
{
  if(is.data.frame(dat))
  {
    stopifnot("SNP" %in% names(dat))
    stopifnot(grepl("rs\\d+",dat$SNP[1]))
    query_snp_list <- dat$SNP
  }
  else
  {
    stopifnot(is.character(dat)&grepl("rs\\d+",dat[1]))
    query_snp_list <- dat
  }
  
  # maximum queries every time is 100
  requery_limit <- ceiling(0.01*length(query_snp_list))
  requery_res <- lapply(1:requery_limit,function(i)
  {
    lower_limit_i <- 100*i - 99
    upper_limit_i <- 100*i
    query_snp_list_i <- query_snp_list[lower_limit_i:upper_limit_i]
    query_snp_list_i <- query_snp_list_i[!is.na(query_snp_list_i)]
    phenoscanner::phenoscanner(snpquery = query_snp_list_i, pvalue = 5e-08)$results
  })
  
  # return a data.frame object
  return(do.call(rbind, requery_res))
}
