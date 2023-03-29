#################################################################################
# I delete the original colocalization template, and upload the most rencet one #
#################################################################################

coloc_test <- function(exposure_dat,
                       outcome_dat,
                       type_exposure = "quant",
                       col_pvalues_exposure = "pval",
                       col_N_exposure = "samplesize",
                       col_MAF_exposure = "maf",
                       col_beta_exposure = "beta",
                       col_se_exposure = "se",
                       col_snp_exposure = "SNP",
                       sd_exposure = NA,
                       type_outcome = "cc",
                       col_pvalues_outcome = "pval.outcome",
                       col_N_outcome = "samplesize.outcome",
                       col_MAF_outcome = NA,
                       col_beta_outcome = "beta.outcome",
                       col_se_outcome = "se.outcome",
                       col_snp_outcome = "SNP",
                       prevalence_outcome = NA)
{
  cols_exposure <- c(col_pvalues_exposure,
                     col_N_exposure,
                     col_MAF_exposure,
                     col_beta_exposure,
                     col_se_exposure,
                     col_snp_exposure)
  cols_exposure <- cols_exposure[!is.na(cols_exposure)]
  
  cols_outcome <- c(col_pvalues_outcome,
                    col_N_outcome,
                    col_MAF_outcome,
                    col_beta_outcome,
                    col_se_outcome,
                    col_snp_outcome)
  cols_outcome <- cols_outcome[!is.na(cols_outcome)]
  
  stopifnot(all(cols_exposure %in% names(exposure_dat)))
  stopifnot(all(cols_outcome %in% names(outcome_dat)))
  
  snp_overlap <- intersect(exposure_dat[[col_snp_exposure]],
                           outcome_dat[[col_snp_outcome]])
  snp_overlap <- unique(snp_overlap)
  
  exposure_dat <- exposure_dat[exposure_dat[[col_snp_exposure]] %in% snp_overlap,]
  outcome_dat <- outcome_dat[outcome_dat[[col_snp_outcome]] %in% snp_overlap,]
  
  exposure_dat <- exposure_dat[order(exposure_dat[[col_snp_exposure]]),]
  outcome_dat <- outcome_dat[order(outcome_dat[[col_snp_outcome]]),]
  
  exposure_list <- list()
  outcome_list <- list()
  
  for (i in 1:8) {
    list_element <- c("pvalues","N","MAF",
                      "beta","varbeta","snps",
                      "type","sdY")[i]
    col_element <- c(col_pvalues_exposure, col_N_exposure, col_MAF_exposure,
                     col_beta_exposure, col_se_exposure, col_snp_exposure,
                     type_exposure, sd_exposure)[i]
    if(!is.na(col_element)){
      if(list_element %in% c("type","sdY"))
      {
        if(list_element=="sdY") col_element <- as.numeric(col_element)
        exposure_list[[list_element]] <- col_element
      }
      else
      {
        if(list_element == "varbeta")
        {
          exposure_list[[list_element]] <- exposure_dat[[col_element]]*exposure_dat[[col_element]]
        }
        else
        {
          exposure_list[[list_element]] <- exposure_dat[[col_element]]
        }
      }
    }
  }
  
  for (j in 1:8) {
    list_element <- c("pvalues","N","MAF",
                      "beta","varbeta","snps",
                      "type","s")[j]
    col_element <- c(col_pvalues_outcome, col_N_outcome, col_MAF_outcome,
                     col_beta_outcome, col_se_outcome, col_snp_outcome,
                     type_outcome, prevalence_outcome)[j]
    if(!is.na(col_element)){
      if(list_element %in% c("type","s"))
      {
        if(list_element=="s") col_element <- as.numeric(col_element)
        outcome_list[[list_element]] <- col_element
      }
      else
      {
        if(list_element == "varbeta")
        {
          outcome_list[[list_element]] <- outcome_dat[[col_element]]*outcome_dat[[col_element]]
        }
        else
        {
          outcome_list[[list_element]] <- outcome_dat[[col_element]]
        }
      }
    }
  }
  
  coloc::coloc.abf(exposure_list,outcome_list, p1 = 1e-04, p2 = 1e-04,p12 = 1e-05)
}
