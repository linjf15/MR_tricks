# Post-GWAS analysis and Mendelian Randomization
During the study based on Mendelian Randomization, several functions were proposed to solve several common problems, including: 
* how to get rsID from chromosome and position: `get_rsID_from_chr_pos`;
* how to find allele frequency, proxy of a given SNP: `snp_add_eaf`, `find_proxy`;
* how to quickly conduct Phenoscanner for a given list of SNP: `snp_phenoscanner`;
* how to harmonize exposure data and local outcome data, if direct removal is not allowed and proxy is required: `snp_replace_proxy`, `harmonise_data_local`
* more functions are undering building
