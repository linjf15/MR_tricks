mr_meta <- function(dat_list,
                    mr_method = c("Wald ratio","Inverse variance weighted"),
                    meta_method = "FE",
                    output = "meta")
{
  requireNamespace("metafor", quietly = TRUE)
  stopifnot(output %in% c("meta","summary"))
  
  mr_method_list <- get_name_from_abbr(mr_method = mr_method)
  
  dat_method_filtering <- lapply(1:length(dat_list),
    function(i) {
      dat_with_method <- dat_list[[i]][dat_list[[i]]$method %in% mr_method_list,]
      dat_with_method$cohort <- i
      dat_with_method
      })
  dat_method_filtering <- do.call(rbind, dat_method_filtering)
  
  exposure_list <- unique(dat_method_filtering$exposure)
  
  if(output=="meta")
  {
    res_meta <- lapply(1:length(exposure_list),function(i)
    {
      dat_i <- dat_method_filtering[dat_method_filtering$exposure %in% exposure_list[i],]
      res_i <- metafor::rma(yi = dat_i$b,sei = dat_i$se, method = meta_method)
      data.frame(exposure = exposure_list[i],
                 b = res_i$beta[1],
                 se = res_i$se,
                 pval = res_i$pval,
                 meta = nrow(dat_i),
                 I2 = res_i$I2,
                 p_het = res_i$QEp)
    })
  }
  if(output=="summary")
  {
    res_meta <- lapply(1:length(exposure_list),function(i)
    {
      dat_i <- dat_method_filtering[dat_method_filtering$exposure %in% exposure_list[i],]
      res_i <- metafor::rma(yi = dat_i$b,sei = dat_i$se, method = meta_method)
      res_i <- data.frame(exposure = exposure_list[i],
                          b = res_i$beta[1],
                          se = res_i$se,
                          pval = res_i$pval,
                          cohort = "meta")
      rbind(dat_i[,c("exposure","b","se","pval","cohort")],
            res_i)
    })
  }
  
  res_meta <- do.call(rbind,res_meta)
  res_meta$lo_ci <- res_meta$b - 1.96 * res_meta$se
  res_meta$up_ci <- res_meta$b + 1.96 * res_meta$se
  res_meta$or <- exp(res_meta$b)
  res_meta$or_lci95 <- exp(res_meta$lo_ci)
  res_meta$or_uci95 <- exp(res_meta$up_ci)
  return(res_meta)
}
