# The original version should do some manually working
# With more understandings about dbsnp (https://www.ncbi.nlm.nih.gov/snp/) and ensembl (http://rest.ensembl.org/), I re-wrote the function
# The function provides two options for database, one for dbsnp, another for ensembl. I prefer dbsnp and add additional function.
# Also, build is considered

get_rsID_from_chrpos <- function(dat = Yang2021_CSF_Marker[31:40,],
                                col_chr = "chr",
                                col_pos = "start",
                                col_ref_allele = "refAllele",
                                col_alt_allele = "altAllele",
                                build = "37",
                                database = "dbsnp")
{
  library(tidyverse)
  
  stopifnot(database %in% c("ensembl","dbsnp"))
  stopifnot(build %in% c("37","38"))
  
  if(database == "ensembl")
  {
    # Create and get a url
    server <- ifelse(build == "37",
                     "http://grch37.rest.ensembl.org",
                     "http://rest.ensembl.org")
    
    query_term <- paste0(
      server,"/vep/human/region/",dat[[col_chr]],":",
      dat[[col_pos]],"-",dat[[col_pos]],"/",dat[[col_alt_allele]],"?")
    
    query_term_alt <- paste0(
      server,"/vep/human/region/",dat[[col_chr]],":",
      dat[[col_pos]],"-",dat[[col_pos]],"/",dat[[col_ref_allele]],"?")
    
    SNP <- lapply(1:nrow(dat), function(i)
    {
      print(paste0("searching for No. ", i, " SNP"))
      query_res <- httr::GET(query_term[i], 
                             httr::content_type("application/json"))
      
      httr::warn_for_status(query_res)
      
      # Convert R objects from JSON
      query_res <- httr::content(query_res)
      res_df <- jsonlite::fromJSON(jsonlite::toJSON(query_res))
      snp <- res_df$colocated_variants[[1]][["id"]][[1]]
      if(is.null(snp))
      {
        query_res <- httr::GET(query_term_alt[i], 
                               httr::content_type("application/json"))
        
        httr::warn_for_status(query_res)
        
        # Convert R objects from JSON
        query_res <- httr::content(query_res)
        res_df <- jsonlite::fromJSON(jsonlite::toJSON(query_res))
        snp <- res_df$colocated_variants[[1]][["id"]][[1]]
        if(is.null(snp)) return(NA)
        else return(snp)
      }
      else return(snp)
    }
    )
    
    dat$SNP <- unlist(SNP)
  }
  
  
  if(database == "dbsnp")
  {
    search_build <- ifelse(build == "37","[POSITION_GRCH37]","[Base Position]")
    
    query_term <- paste0(dat[[col_chr]],"[CHR] AND Homo[ORGN] AND ",
                         dat[[col_pos]],search_build)
    
    SNP <- lapply(1:nrow(dat), function(i)
    {
      print(paste0("searching for No. ", i, " SNP"))
      snp <- unlist(rentrez::entrez_search(db="snp", term=query_term[i])$ids)
      if(is.null(snp)) return(NA)
      else return(paste0("rs",snp[length(snp)]))
    }
    )
    
    dat$SNP <- unlist(SNP)
    
    if(nrow(dat) != length(SNP[!is.na(SNP)]))
    {
      dat_with_snp <- dat[!is.na(dat$SNP),]
      dat_without_snp <- dat[is.na(dat$SNP),]
      query_term_alt <- paste0(dat_without_snp[[col_chr]],"[CHR] AND Homo[ORGN] AND ",
                               dat_without_snp[[col_pos]]+1,search_build)
      SNP <- lapply(1:nrow(dat_without_snp), function(i)
      {
        print(paste0("Researching for No. ", i, " SNP"))
        snp <- unlist(rentrez::entrez_search(db="snp", term=query_term_alt[i])$ids)
        if(is.null(snp)) return(NA)
        else return(paste0("rs",snp[length(snp)]))
      }
      )
      dat_without_snp$SNP <- unlist(SNP)
      dat <- rbind(dat_with_snp, dat_without_snp)
    }
  }
  
  dat
}
