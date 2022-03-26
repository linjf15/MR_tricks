# Post-GWAS analysis and Mendelian Randomization
During the study based on Mendelian Randomization, several functions were proposed to solve several common problems, including: 

## Preprocessing for exposure GWASs
* how to get rsID from chromosome and position: `get_rsID_from_chr_pos`;
* how to find allele frequency, proxy of a given SNP: `snp_add_eaf`, `snp_add_eaf_manual`,`find_proxy`;
* how to search pleiotropy for a given SNP: `snp_phenoscanner`;

## Exposure and outcome data harmonization
* how to harmonize if direct removal is not allowed and proxy is required: `snp_replace_proxy`, `harmonise_data_modified`

## Main analysis for Mendelian Randomization
* how to conduct meta-analysis to pool estimates from different cohorts: `mr_meta`
* `get_cohort_from_number` transforms a number into a name of cohort
* `get_name_from_abbr` transforms an abbreviation version of MR method into a full name

## More functions are under developing...
