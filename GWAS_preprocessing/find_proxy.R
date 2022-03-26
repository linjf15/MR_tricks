find_proxy <- function(snp = "rs10001", r2_threshold = 0.8, build = "37", pop = "EUR")
{
  stopifnot(build %in% c("37","38"))
  
  # Create and get a url
  server <- ifelse(build == "37","http://grch37.rest.ensembl.org","http://rest.ensembl.org")
  ext <- paste0("/ld/human/",snp,"/1000GENOMES:phase_3:", pop)
  url <- paste(server, ext, sep = "")
  print(paste0("searching proxies for ",snp," ......"))
  res <- httr::GET(url, httr::content_type("application/json"))
  
  # Converts http errors to R errors or warnings
  httr::stop_for_status(res)
  
  # Convert R objects from JSON
  res <- httr::content(res)
  proxy_df <- jsonlite::fromJSON(jsonlite::toJSON(res))
  proxy_df <- lapply(proxy_df, unlist)
  
  # Extract SNP associated with targeted SNP with r2 > r2_threshold
  snp_r2 <- proxy_df$r2[proxy_df$r2>r2_threshold]
  snp_v <- proxy_df$variation2[proxy_df$r2>r2_threshold]
  
  # If the query returns nothing, the output would be NULL
  # Else, returning a ordered vectors of proxy SNPs
  if(is.null(snp_v)) return(NULL)
  else return(snp_v[order(snp_r2, decreasing=T)])
}
