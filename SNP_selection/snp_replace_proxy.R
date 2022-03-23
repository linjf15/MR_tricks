snp_replace_proxy <- function(dat, snp_proxy, type = "exposure", build = "37", pop = "EUR")
{
  stopifnot(nrow(dat) == 1)
  stopifnot("SNP" %in% names(dat))
  stopifnot(paste0("effect_allele.",type) %in% names(dat))
  stopifnot(paste0("eaf.",type) %in% names(dat))
  stopifnot(build %in% c("37","38"))
  
  # Create and get a url
  server <- ifelse(build == "37","http://grch37.rest.ensembl.org","http://rest.ensembl.org")
  pop <- paste0("1000GENOMES:phase_3:",pop)
  
  ext <- paste0("/variation/Homo_sapiens/", snp_proxy, "?content-type=application/json;pops=1")
  url <- paste(server, ext, sep = "")
  res <- httr::GET(url)
  
  # Converts http errors to R errors or warnings
  httr::stop_for_status(res)
  
  # Convert R objects from JSON
  res <- httr::content(res)
  res_pop <- jsonlite::fromJSON(jsonlite::toJSON(res))$populations
  
  # Filter query results based on population set
  res_pop <- res_pop[res_pop$population == pop,]
  
  ref_allele <- res_pop[1,"allele"][[1]]
  alt_allele <- res_pop[2,"allele"][[1]]
  ref_af <- res_pop[1,"frequency"][[1]]
  
index <- ifelse(
    (dat[[paste0("eaf.",type)]]<0.5&ref_af<0.5)|(dat[[paste0("eaf.",type)]]>0.5&ref_af>0.5),
    1,2
  )

  dat[["SNP"]] <- snp_proxy
  dat[[paste0("effect_allele.",type)]] <- c(ref_allele, alt_allele)[index]
  dat[[paste0("eaf.",type)]] <- abs(1 - index + ref_af)
    
  if(paste0("other_allele.",type) %in% names(dat))
  {
    dat[[paste0("other_allele.",type)]] <- c(ref_allele, alt_allele)[3 - index]
  }

  return(dat)
}
