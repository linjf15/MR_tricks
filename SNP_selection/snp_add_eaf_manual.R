# snp_add_eaf generally provides a raw path to extract effect allele frequency
# However, some DEL/IN SNP should be further checked

snp_add_eaf_manual <- function(dat, 
                               old_SNP = "rs1484845180", new_SNP = "rs201903842", 
                               force_replace = FALSE,
                               build = "37", pop = "EUR")
{
  requireNamespace("stringr",quietly = TRUE)
  stopifnot(build %in% c("37","38"))
  stopifnot("SNP" %in% names(dat))
  dat_i <- dat[dat$SNP == old_SNP,]
  stopifnot(nrow(dat_i)==1)
  dat_i$SNP <- new_SNP
  
  # Create and get a url
  server <- ifelse(build == "37","http://grch37.rest.ensembl.org","http://rest.ensembl.org")
  pop <- paste0("1000GENOMES:phase_3:",pop)
  
  ext <- paste0("/variation/Homo_sapiens/",new_SNP, "?content-type=application/json;pops=1")
  url <- paste(server, ext, sep = "")
  res <- httr::GET(url)
  
  # Converts http errors to R errors or warnings
  httr::stop_for_status(res)
  
  # Convert R objects from JSON
  res <- httr::content(res)
  res_pop <- jsonlite::fromJSON(jsonlite::toJSON(res))$populations
  
  # Filter query results based on population set
  res_pop <- try(res_pop[res_pop$population == pop,])
  if("try-error" %in% class(res_pop))
  {
    print(paste0("There is not information for population ",pop))
    queried_effect_allele <- "NR"
    queried_other_allele <- "NR"
    queried_eaf <- -1
  }
  else
  {
    if(nrow(res_pop)==0)
    {
      print(paste0("There is not information for population ",pop))
      queried_effect_allele <- "NR"
      queried_other_allele <- "NR"
      queried_eaf <- -1
    }
    else
    {
      queried_effect_allele <- res_pop[1,"allele"][[1]]
      queried_other_allele <- res_pop[2,"allele"][[1]]
      queried_eaf <- res_pop[1,"frequency"][[1]]    
    }
  }
  
  
  effect_allele <- ifelse("effect_allele.exposure" %in% names(dat),
                          dat_i$effect_allele.exposure,
                          dat_i$effect_allele)
  
  other_allele <- ifelse("effect_allele.exposure" %in% names(dat),
                         dat_i$other_allele.exposure,
                         dat_i$other_allele)
  
  if(force_replace)
  {
    len_effect_allele <- nchar(effect_allele)
    len_other_allele <- nchar(other_allele)
    len_queried_effect_allele <- nchar(queried_effect_allele)
    len_queried_other_allele <- nchar(queried_other_allele)
    
    if((len_queried_effect_allele-len_effect_allele) == 
       (len_queried_other_allele-len_other_allele))
    {
      dat_i$eaf.exposure <- queried_eaf
      dat_i$eaf <- dat_i$eaf.exposure 
      dat_i$reliability.exposure <- "high"
    }
    if((len_queried_effect_allele-len_other_allele) == 
       (len_queried_other_allele-len_effect_allele))
    {
      dat_i$eaf.exposure <- 1 - queried_eaf
      dat_i$eaf <- dat_i$eaf.exposure 
      dat_i$reliability.exposure <- "high"
    }
    
  }
  else
  {
    if(str_detect(effect_allele,queried_effect_allele))
    {
      print(paste0(queried_effect_allele, " in ",effect_allele))
      dat_i$eaf.exposure <- queried_eaf
      dat_i$eaf <- dat_i$eaf.exposure 
      dat_i$reliability.exposure <- "high"
    }
    if(str_detect(other_allele,queried_effect_allele))
    {
      print(paste0(queried_effect_allele, " in ",other_allele))
      dat_i$eaf.exposure <- 1 - queried_eaf
      dat_i$eaf <- dat_i$eaf.exposure 
      dat_i$reliability.exposure <- "high"
    }
  }
  return(rbind(dat_i[names(dat)],
               dat[dat$SNP != old_SNP,][names(dat)]))
}
