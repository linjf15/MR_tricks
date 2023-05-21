# Post-GWAS analysis and Mendelian Randomization
During the study based on Mendelian Randomization, several functions were proposed to solve several common problems, including: 

## Special note
Recently, I noticed that some uploaders in Bilibili seem to sell my codes, such as [Shining94](https://space.bilibili.com/295917932/). The original purpose of this project is to share my solutions to certain problems in Mendelian Randomization without charge. I urge those guys stop selling my codes, and I decide to not share my recent findings on MR. If you have any questions, leave me a message.

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
* [MR dictionary](https://mr-dictionary.mrcieu.ac.uk/)
* Cohorts: [FinnGen](https://finngen.gitbook.io/documentation/), [UK Biobank](https://www.leelabsg.org/resources)
* Algorithm: [MR-RAPS](https://github.com/qingyuanzhao/mr.raps), [MR-PRESSO](https://github.com/rondolab/MR-PRESSO), [MR-Mix](https://github.com/gqi/MRMix), [IMRP](https://github.com/xiaofengzhucase/IMRP), [MR-GENIUS](https://github.com/bluosun/MR-GENIUS), [MR-AID](https://github.com/yuanzhongshang/MRAID), [MR-TRYX](https://github.com/explodecomputer/tryx), [NLMR](https://github.com/jrs95/nlmr), [MVMR](https://github.com/WSpiller/MVMR), [BayesMR](https://github.com/igbucur/BayesMR), [CFMR](https://github.com/william-denault/CFMR), [OMR](https://github.com/wanglu205/OMR), [MR-JTI](https://github.com/gamazonlab/MR-JTI)
* Packages: [ieugwasr](https://mrcieu.github.io/ieugwasr/articles/)
