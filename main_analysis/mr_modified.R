mr_modified <- function (dat, 
                         parameters = default_parameters(), 
                         method_list = subset(mr_method_list(), use_by_default)$obj) 
{
  calculate_prop_var_explained <- function(dat)
  {
    prop_var_explained_sum <- 0
    for (i in 1:nrow(dat)) {
      dat_i <- dat %>% dplyr::filter(row_number()==i)
      study_id <- dat_i$gene.exposure %>% str_split("_") %>% map(1) %>% unlist()
      sample_size <- case_when(study_id=="Sun2018"~3301,
                               study_id=="Emilsson2018"~3200,
                               study_id=="Suhre2017"~1000,
                               study_id=="Yao2018"~6861,
                               study_id=="Folkersen2020"~21758,
                               study_id=="Yang2021"~835)
      prop_var_explained_i <- dat_i %>% 
        dplyr::mutate(prop_var_explained = (beta.exposure*beta.exposure)/(beta.exposure*beta.exposure+se.exposure*se.exposure*sample_size)) %>% 
        dplyr::pull(prop_var_explained)
      prop_var_explained_sum <- prop_var_explained_sum + prop_var_explained_i
    }
    
    return(prop_var_explained_sum)
  }
  mr_raps_modified <- function (b_exp, b_out, se_exp, se_out,parameters) 
  {
    out <- try(suppressMessages(mr.raps::mr.raps(b_exp, b_out, se_exp, se_out,
                                                 over.dispersion = parameters$over.dispersion, 
                                                 loss.function = parameters$loss.function,
                                                 diagnosis = FALSE)),
               silent = T)
    if ('try-error' %in% class(out))
    {
      output = list(b = NA, se = NA, pval = NA, nsnp = NA)
    }
    else
    {
      output = list(b = out$beta.hat, se = out$beta.se, pval = pnorm(-abs(out$beta.hat/out$beta.se)) * 
                      2, nsnp = length(b_exp))
    }
    return(output)
  }
  method_list_modified <- method_list %>%
    str_replace_all("mr_raps","mr_raps_modified")
  
  mr_tab <- plyr::ddply(dat, c("id.exposure", "id.outcome"), 
                        function(x1) {
                          x <- subset(x1, mr_keep)
                          
                          if (nrow(x) == 0) {
                            message("No SNPs available for MR analysis of '", 
                                    x1$id.exposure[1], "' on '", x1$id.outcome[1], 
                                    "'")
                            return(NULL)
                          }
                          else {
                            message("Analysing '", x1$id.exposure[1], 
                                    "' on '", x1$id.outcome[1], "'")
                          }
                          res <- lapply(method_list_modified, function(meth) {
                            get(meth)(x$beta.exposure, x$beta.outcome, 
                                      x$se.exposure, x$se.outcome, 
                                      parameters)
                          })
                          
                          methl <- mr_method_list()
                          mr_tab <- data.frame(outcome = x$outcome[1], 
                                                exposure = x$exposure[1], 
                                                method = methl$name[match(method_list, methl$obj)], 
                                                nsnp = sapply(res, function(x) x$nsnp), 
                                                b = sapply(res, function(x) x$b), 
                                                se = sapply(res, function(x) x$se), 
                                                pval = sapply(res, function(x) x$pval))
                          
                          mr_tab <- subset(mr_tab, !(is.na(b) & is.na(se) & 
                                                        is.na(pval)))
                          
                          # proportion of variance explained
                          
                          prop_var_explained <- try(calculate_prop_var_explained(x))
                          if ('try-error' %in% class(prop_var_explained)) prop_var_explained <- NA
                          mr_tab$pve <- prop_var_explained
                          
                          return(mr_tab)
                        })
  return(mr_tab)
}
