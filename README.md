# Post-GWAS analysis and Mendelian Randomization
During the study based on Mendelian Randomization, several functions were proposed to solve several common problems, including: 

## Preprocessing for exposure GWASs
* how to get rsID from chromosome and position: `get_rsID_from_chr_pos`
* how to find allele frequency, proxy of a given SNP: `snp_add_eaf`, `snp_add_eaf_manual`,`find_proxy`
* how to search pleiotropy for a given SNP: `snp_phenoscanner`

## Exposure and outcome data harmonization
* how to harmonize if direct removal is not allowed and proxy is required: `snp_replace_proxy`, `harmonise_data_modified`

## Main analysis for Mendelian Randomization
* `mr_meta` is helpful conduct meta-analysis to pool estimates from different cohorts
* `mr_modified` is alternative function to `mr` in TwoSampleMR package, since the mr-raps cannot not be conducted in the original function
* `get_cohort_from_number` transforms a number into a name of cohort
* `get_method_from_abbr` transforms an abbreviation version of MR method into a full name
* `get_protein_from_uniprot` transforms a UniProt ID to full name of protein, with UniProtRefPanel.xlsx as the reference panel

## More functions are under developing...
* data visualization: volcano plot, forest plot...
* colocolization
* TWAS
* ...

## Useful references
* [Usage of ieugwasr](https://mrcieu.github.io/ieugwasr/articles/)
* [FinnGen](https://finngen.gitbook.io/documentation/)
* [UK Biobank](https://www.leelabsg.org/resources)
* [MR dictionary](https://mr-dictionary.mrcieu.ac.uk/)
