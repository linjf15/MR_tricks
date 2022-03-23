snp_add_eaf <- function(dat, build = "37", pop = "EUR")
{
  stopifnot(build %in% c("37","38"))
  stopifnot("SNP" %in% names(exposure_MAP))
  
  # Create and get a url
  server <- ifelse(build == "37","http://grch37.rest.ensembl.org","http://rest.ensembl.org")
  pop <- paste0("1000GENOMES:phase_3:",pop)
  
  res_tab <- lapply(1:nrow(dat), function(i)
  {
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
    res_pop <- res_pop[res_pop$population == pop,]
    allele <- res_pop[1,"allele"][[1]]
    af <- res_pop[1,"frequency"][[1]]
    
    if("effect_allele.exposure" %in% names(dat))
    {
      dat_i$eaf.exposure <- ifelse(dat_i$effect_allele.exposure == allele,af,1-af)
    }
    if("effect_allele" %in% names(dat))
    {
      dat_i$eaf <- ifelse(dat_i$effect_allele == allele,af,1-af)
    }
    dat_i[names(dat)]
  })
  
  return(do.call(rbind, res_tab))
}
