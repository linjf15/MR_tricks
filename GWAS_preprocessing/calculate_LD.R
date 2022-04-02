calculate_LD <- function(snp1 = "rs6792369", snp2 = "rs1042779", build = "37", pop = "EUR")
{
  stopifnot(build %in% c("37","38"))

  # Create and get a url
  server <- ifelse(build == "37","http://grch37.rest.ensembl.org","http://rest.ensembl.org")
  ext <- paste0("/ld/human/pairwise/",snp1,"/",snp2,"?population_name=1000GENOMES:phase_3:", pop)
  url <- paste(server, ext, sep = "")
  print(paste0("Calculating LD for ",snp1," and ",snp2," ......"))
  res <- httr::GET(url, httr::content_type("application/json"))
  
  # Converts http errors to R errors or warnings
  httr::stop_for_status(res)
  
  # Convert R objects from JSON
  res <- httr::content(res)
  ld_df <- jsonlite::fromJSON(jsonlite::toJSON(res))
  
  if("r2" %in% names(ld_df)) return(unlist(ld_df$r2))
  else return(url)
}
