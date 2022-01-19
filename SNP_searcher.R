library(tidyverse)

# 1. get chr (CHROM), pos (POS), ref (REF), and alt(ALT) and generate formatted dat

Yang2021_CSF_Marker <- readxl::read_excel("Yang2021_CSF_Marker.xlsx") %>% 
  dplyr::mutate(snpID = str_split(protein_snpID_roundID,"_") %>% 
                  map(2) %>% unlist()) %>%
  dplyr::pull(snpID) %>% 
  unique() %>% 
  tibble() %>% 
  rename(snpID = 1) %>% 
  dplyr::mutate(CHROM = str_split(snpID,":") %>% map(1) %>% unlist(),
                POS = str_split(snpID,":") %>% map(2) %>% unlist(),
                REF = str_split(snpID,":") %>% map(3) %>% unlist(),
                ALT = str_split(snpID,":") %>% map(4) %>% unlist()) %>% 
  dplyr::mutate(CHROM = as.numeric(CHROM),
                POS = as.numeric(POS)) %>% 
  dplyr::select(CHROM,POS,REF,ALT)


# 2.Use CHROM, POS, REF and ALT to generate VCF files
# It is noteworthy that, VEP (variant effect predictor) is sensitive to the first allele given
# If an ALT was present in the REF column, a NULL result may be returned
# We need to reverse the position of REF and ALT

export_VCF <- function(formatted_dat = Yang2021_CSF_Marker)
{
  output_file <- paste0(as.character(substitute(formatted_dat)),"_SNPlist.tsv")
  
  reverse_formatted_dat <- formatted_dat %>% 
    dplyr::mutate(REF_temp = ALT,
                  ALT_temp = REF) %>% 
    dplyr::select(CHROM,POS,REF_temp,ALT_temp) %>% 
    dplyr::rename(REF = REF_temp,
                  ALT = ALT_temp)
  
  write.table(formatted_dat %>% 
                dplyr::bind_rows(reverse_formatted_dat) %>% 
                dplyr::mutate(ID = ".",
                              QUAL = ".",
                              FILTER = ".",
                              INFO = ".") %>% 
                dplyr::select(CHROM,POS,ID,REF,
                              ALT,QUAL,FILTER,INFO) %>% 
                # it is very important to arrange, or VEP would fail
                dplyr::arrange(CHROM,POS),
              output_file,
              sep = " ",
              row.names = F,
              col.names = F,
              quote = F)
}

export_VCF()


# 3. For VCF file, uploaded to http://grch37.ensembl.org/Tools/VEP
# Create a New job, and paste the data in VCF file
# After the job is done, download the result in .txt format and renamed it as "Yang2021_CSF_Marker_supp.txt"


import_VEP <- function(VEP_result = "Yang2021_CSF_Marker_supp.txt")
{
  read_delim(VEP_result, 
             delim = "\t", 
             escape_double = FALSE, 
             trim_ws = TRUE) %>% 
    dplyr::select(Location, Allele,Existing_variation,AF) %>% 
    unique() %>% 
    dplyr::mutate(CHROM = Location %>% str_split(":") %>% map(1) %>% unlist() %>% as.numeric(),
                  POS_START = Location %>% str_split("-") %>% map(1) %>% unlist() %>% 
                    str_split(":") %>% map(2) %>% unlist() %>% as.numeric(),
                  POS_END = Location %>% str_split("-") %>% map(2) %>% unlist() %>% as.numeric(),
                  SNP = Existing_variation %>% str_extract("rs\\d+") %>% unlist()) %>% 
    dplyr::select(-Existing_variation) %>% 
    dplyr::rename(ALLELE=Allele) %>% 
    dplyr::filter(!is.na(SNP)) %>% 
    dplyr::group_by(SNP) %>% 
    dplyr::filter(row_number()==1) %>% 
    dplyr::ungroup()
}

Yang2021_CSF_Marker_supp <- import_VEP()



# criteria 1: A SNP should show up for only once
# create a function to determine if both allele is 
# An allele = "-" may also be of significance, eg. 1:2513652-2513656
# criteria 2: 



add_snp <- function(formatted_dat = Yang2021_CSF_Marker,
                    VEP_dat = Yang2021_CSF_Marker_supp,
                    VEP_output = FALSE,
                    Search = TRUE)
{
  # first, we will classify the data into 2 tiers
  # tier 1: both alleles are in A/T/C/G, and the result is of high accuracy
  # tier 2: at least one allele is indel, and the result might be of low accuracy
  # if needed, you might rewrite this part to ensure output quality
  
  dat_origin <- formatted_dat %>% 
    dplyr::mutate(tier=str_length(REF)+str_length(ALT)) %>% 
    dplyr::mutate(tier=case_when(tier==2~1,
                                 tier>2~2))
  
  # the reference SNPs are in the previous VEP files
  dat_ref <- VEP_dat %>% 
    dplyr::filter(str_length(ALLELE)==1) %>% # this is a strong filter
    dplyr::group_by(CHROM,POS_START) %>% 
    dplyr::arrange(SNP) %>% 
    dplyr::filter(row_number()==1) %>% 
    dplyr::ungroup()
  
  # the key to link two dataset is the chr:pos
  dat_tier_1 <- dat_origin %>% 
    dplyr::filter(tier == 1) %>% 
    dplyr::left_join(dat_ref,
                     by = c("CHROM"="CHROM","POS"="POS_START")
                     )
  
  dat_tier_1_finished <- dat_tier_1 %>% 
    dplyr::filter(!is.na(SNP))
  
  # delete SNP that has been used/ matched and updated the dat_ref
  dat_ref <- VEP_dat %>% 
    dplyr::filter(!SNP %in% dat_tier_1_finished$SNP) %>% # this is a weak filter
    dplyr::group_by(CHROM,POS_START) %>% 
    dplyr::arrange(SNP) %>% 
    dplyr::filter(row_number()==1) %>% 
    dplyr::ungroup()
  
  # the left data was classified as tier_2 dataset
  dat_tier_2 <- dat_tier_1 %>% 
    dplyr::filter(is.na(SNP)) %>%
    dplyr::select(all_of(names(dat_origin))) %>% 
    dplyr::bind_rows(dat_origin %>% filter(tier == 2)) %>% 
    dplyr::left_join(dat_ref,
                     by = c("CHROM"="CHROM","POS"="POS_START")
                     )
  
  dat_tier_2_finished <- dat_tier_2 %>% 
    dplyr::filter(!is.na(SNP))
  
  dat_found <- dat_tier_1_finished %>% 
    bind_rows(dat_tier_2_finished)
  
  dat_not_found <- dat_tier_2 %>% 
    dplyr::filter(is.na(SNP))
  
  if(VEP_output)
  {
    output_file <- paste0(as.character(substitute(formatted_dat)),"_SNPlist2.tsv")
    
    # we would like to generate a second VCF data
    # to confirm that the SNP info was truely missing
    write.table(dat_not_found %>% 
                  dplyr::mutate(ID = ".",
                                QUAL = ".",
                                FILTER = ".",
                                INFO = ".") %>% 
                  dplyr::select(CHROM,POS,ID,REF,
                                ALT,QUAL,FILTER,INFO) %>% 
                  dplyr::arrange(CHROM,POS),
                output_file,
                sep = " ",
                row.names = F,
                col.names = F,
                quote = F)
  }
  
  if(Search)
  {
    dat_not_found <- dat_not_found %>% 
      dplyr::arrange(CHROM,POS) %>% 
      dplyr::mutate(search_term = paste0(CHROM,"[CHR] AND Homo[ORGN] AND ",POS,"[Base Position Previous]"))
  
    return(list(dat_found,dat_not_found))
  }
  else
  {
    return(list(dat_found))
  }
}


# In the final step, we will do two jobs

manual_search_dat <- add_snp(Search = T)[[2]] %>% 
  dplyr::filter(str_length(REF)>=2|str_length(ALT)>=2)

automated_search_dat <- add_snp(Search = T)[[2]] %>% 
  dplyr::filter(str_length(REF)<2&str_length(ALT)<2)

# 1. For REF or ALT not in A/T/C/G, a manual search in https://www.ncbi.nlm.nih.gov/snp/ is highly recommended

# 19 55064008 rs1380947605
# 19 55132712 rs34781265

# 2. For both REF and ALT in A/T/C/G, a automated search could be conducted by rentrez


library(rentrez)

search_snp <- function(search_dat = add_snp()[[2]])
{
  dat_supp <- tibble()
  for (each_term in search_dat$search_term) {
    snp_search <- entrez_search(db="snp", 
                                term=each_term)
    snp_search_res <- snp_search$ids
    if(length(snp_search_res)==0) snp_id = NA
    else
    {
      snp_id <- paste0("rs",snp_search_res[length(snp_search_res)])
    }
    dat_supp <- dat_supp %>% 
      bind_rows(c(search_term = each_term, SNP = snp_id))
    # print(snp_id)
    Sys.sleep(0.1)
  }
  
  dat_supp %>% 
    dplyr::bind_cols(search_dat) %>% 
    dplyr::select(SNP,CHROM,POS)
}


Yang2021_CSF_Marker_entrez <- search_snp(search_dat = automated_search_dat)

Yang2021_CSF_Marker_Final <- automated_search_dat %>% 
  dplyr::select(-SNP,-search_term) %>% 
  dplyr::left_join(Yang2021_CSF_Marker_entrez,
                   by = c("CHROM"="CHROM","POS"="POS")) %>% 
  dplyr::mutate(SNP = case_when(CHROM==19&POS==55064008~"rs1380947605",
                                CHROM==19&POS==55132712~"rs34781265",
                                TRUE~SNP)) %>% 
  dplyr::bind_rows(add_snp(Search = F))


# Finally, this work is done

Yang_CSF_Original <- readxl::read_excel("Yang_CSF_Original.xlsx")


writexl::write_xlsx(Yang_CSF_Original %>% 
  dplyr::left_join(Yang2021_CSF_Marker_Final %>% 
                     dplyr::select(CHROM,POS,REF,ALT,SNP),
                   by = c("chr"="CHROM","start" = "POS", "refAllele"="REF","altAllele"="ALT")
                   ),
  "Yang_CSF_Updated.xlsx"
  )

Yang_Blood_Original <- readxl::read_excel("Yang_Blood_Original.xlsx")

writexl::write_xlsx(Yang_Blood_Original %>% 
                      dplyr::left_join(Yang2021_CSF_Marker_Final %>% 
                                         dplyr::select(CHROM,POS,REF,ALT,SNP),
                                       by = c("chr"="CHROM","start" = "POS", "refAllele"="REF","altAllele"="ALT")
                      ),
                    "Yang_Blood_Updated.xlsx"
)

