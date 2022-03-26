### several functions are required for harmonizing data
### TwoSampleMR package
### find_proxy
### snp_replace_proxy

harmonise_data_local <- function(exposure_dat, outcome_dat, r2_thershold = 0.8, build = "37", pop = "EUR")
{
  stopifnot(build %in% c("37","38"))
  
  library(TwoSampleMR)
  # TwoSampleMR::harmonise_data should be used
  # See: https://mrcieu.github.io/TwoSampleMR/articles/harmonise.html
  # I honor the work conducted by Gibran Hemani et al. sincerely

  outcome_from <- ifelse(is.character(outcome_dat),"ieugwas","local")
  outcome_id <- ifelse(is.character(outcome_dat),outcome_dat,NA)
  snp_in_exposure_dat <- exposure_dat$SNP
  if(outcome_from == "ieugwas")
  {
    # palindromes is stated because palindromes would not be simply deleted by the function
    # on the other hand, TwoSampleMR::harmonise_data would comprehensively analyzed it
    outcome_dat <- extract_outcome_data(snp_in_exposure_dat, outcome_id,
                                        proxies = TRUE, rsq = r2_thershold,
                                        palindromes = 1)
  }
  harmonised_dat <- harmonise_data(exposure_dat, outcome_dat)
  harmonised_dat <- harmonised_dat[harmonised_dat$mr_keep == "TRUE",]
  harmonised_dat$data_source.exposure <- "input"
  
  # Get a list of unavailable snps in the outcome GWAS
  all_snp <- snp_in_exposure_dat[grepl("rs\\d+", snp_in_exposure_dat)]
  snp_available <- harmonised_dat$SNP[grepl("rs\\d+", harmonised_dat$SNP)]
  snp_na <- all_snp[!all_snp %in% snp_available]
  
  for (missed_snp in snp_na)
  {
    # For a given unavailable snp, find all the proxy snps with r2 > r2_thershold
    missed_snp_proxy <- find_proxy(missed_snp, r2_thershold, build, pop)
    
    # If there is not proxy for snp unavailable, passed
    if(length(missed_snp_proxy)==0) next
    
    else
    {
      time_search <- 0
      Sys.sleep(0.1)
      print(paste0("Proxies for ",missed_snp, ": ", paste0(missed_snp_proxy,collapse = ", ")))
      
      # Get the exposure dat for the unavailable
      exposure_dat_i <- exposure_dat[exposure_dat$SNP == missed_snp,]

      # Looping for all the proxy snps for the given unavailable snp
      for (proxy_snp_i in missed_snp_proxy)
      {
        # Limiting the times of searching
        time_search <- time_search + 1
        if(time_search > length(missed_snp_proxy)) break
        
        # Searching GWAS result for proxy snp in the outcome data
        if(outcome_from == "ieugwas") 
        {
          outcome_dat_i <- extract_outcome_data(proxy_snp_i,outcome_id,proxies = F)
          if(is.null(outcome_dat_i)) next
        }
        if(outcome_from == "local") 
        {
          if(!proxy_snp_i %in% outcome_dat$SNP) next
          else
          {
            outcome_dat_i <- outcome_dat[outcome_dat$SNP == proxy_snp_i,]
            outcome_dat_i <- snp_replace_proxy(outcome_dat_i, missed_snp, type = "outcome", build, pop)
          }
        }
        
        # Harmonizing data based on the proxy-outcome GWAS
        dat_supp_i <- harmonise_data(exposure_dat_i, outcome_dat_i)
        if(nrow(dat_supp_i)==0) next
        else
        {
          if(dat_supp_i$mr_keep=="TRUE")
          {
            dat_supp_i$data_source.exposure <- proxy_snp_i
            harmonised_dat <- rbind(harmonised_dat, dat_supp_i)
            break
          }
        }
      }
    }
  }
  
  return(harmonised_dat)
}
