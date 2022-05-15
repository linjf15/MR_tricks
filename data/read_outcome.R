######################################################################
# This script documents reading functions for UK Biobank and FinnGen #
######################################################################

library(tidyverse)
library(TwoSampleMR)


read_lee <- function(file_name,file_storage = "D:/R/learning/MR/Outcome_GWAS/UKB_Lee/")
{
  file_path <- paste0(file_storage,file_name,".vcf")
  read_table2(file_path)
}


format_lee <- function(dat,outcome,id)
{
  dat %>% 
    dplyr::select(ID,"#CHROM",POS,ALT,REF,af,beta,sebeta,pval,num_cases,num_controls) %>% 
    dplyr::rename(SNP = ID, 
                  chr = "#CHROM", pos = POS,
                  effect_allele = ALT, other_allele = REF, eaf = af,
                  se = sebeta,
                  ncase = num_cases, ncontrol = num_controls) %>% 
    dplyr::mutate(Phenotype := !!outcome, id := !!id) %>% 
    TwoSampleMR::format_data(type = "outcome",min_pval = 1e-400)
}


read_finn <- function(file_name,file_storage = "D:/R/learning/MR/Outcome_GWAS/FinnGen")
{
  file_path <- paste0(file_storage,"/",file_name,"/",file_name)
  read_delim(file_path,delim = "\t", escape_double = FALSE,trim_ws = TRUE)
}


format_finn <- function(dat,outcome,id,ncase,ncontrol)
{
  dat %>% 
    dplyr::select(rsids,beta,sebeta,af_alt,alt,ref,pval,"#chrom",pos,nearest_genes) %>% 
    dplyr::rename(SNP=rsids,se=sebeta,eaf=af_alt,effect_allele =alt,other_allele=ref,chr = "#chrom",gene = nearest_genes) %>% 
    dplyr::mutate(Phenotype := !!outcome, id := !!id) %>% 
    dplyr::mutate(ncase = ncase, ncontrol = ncontrol) %>% 
    TwoSampleMR::format_data(type = "outcome",min_pval = 1e-400)
}

