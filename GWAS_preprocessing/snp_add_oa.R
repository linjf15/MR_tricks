snp_add_oa <- function(dat, build = "37", pop = "EUR")
{
  stopifnot(build %in% c("37","38"))
  stopifnot("SNP" %in% names(dat))
  
  # Create and get a url
  server <- ifelse(build == "37","http://grch37.rest.ensembl.org","http://rest.ensembl.org")
  pop <- paste0("1000GENOMES:phase_3:",pop)
  name_other_allele <- ifelse("effect_allele.exposure" %in% names(dat),
                              "other_allele.exposure",
                              "other_allele")
  name_effect_allele <- ifelse("effect_allele.exposure" %in% names(dat),
                              "effect_allele.exposure",
                              "effect_allele")
  
  res_tab <- lapply(1:nrow(dat), function(i)
  {
    print(paste0("searching for No.", i, " SNP"))
    dat_i <- dat[i,]

    
    ext <- paste0("/variation/Homo_sapiens/",dat_i$SNP, "?content-type=application/json;pops=1")
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
      dat_i[[name_other_allele]] <- "NR"
      dat_i$ref_effect_allele <- "NR"
    }
    else
    {
      if(nrow(res_pop)!=2)
      {
        print(paste0("There is not information for population ",pop))
        dat_i[[name_other_allele]] <- "NR"
        dat_i$ref_effect_allele <- "NR"
      }
      if(nrow(res_pop)==2)
      {
        queried_effect_allele <- res_pop[1,"allele"][[1]]
        queried_other_allele <- res_pop[2,"allele"][[1]]
        dat_i[[name_other_allele]] <- ifelse(
          dat_i[[name_effect_allele]] == queried_effect_allele,
          queried_other_allele,queried_effect_allele)
        dat_i$ref_effect_allele <- ifelse(
          dat_i[[name_effect_allele]] == queried_effect_allele,
          queried_effect_allele,queried_other_allele)
      }
    }
    dat_i
  })
  return(res_tab)
}
