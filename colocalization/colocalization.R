## This script runs colocalisation analysis between eQTL data for genes near ACE
## with schizophrenia GWAS data

## Input files:
## 1. eQTL data in SMR query format (provided in the smr_query_files folder)
## 2. GWAS data gcta_format
## 3. List of genes for which analysis needs to be done (ENSG_IDs_for_genes_near_ACE.txt)

## Packages required
library(coloc)

##################### 1. EXTRACT EQTL SUMMARY DATA FOR GENES OF INTEREST USING SMR  #######
## NOTE: eqtl summary files are provided in the smr_query_files folder and this step can be skipped
##########################################################################################
# SMR command
system("smr --beqtl-summary myeqtl_data --query 1 --gene gene_id --out myquery")
# where myeqtl_data is the summary-level data from a eQTL study in smr format which can be
# downloaded from https://cnsgenomics.com/software/smr/#DataResource;
# --query sets the p-value threshold. To extract all SNPs use --query 1;
# --gene is the gene_id as provided in the eqtl data. For eQTLGen data this is the ensembl ID
##########################################################################################



##################### 1. COLOC ANALYSIS OF EQTL AND GWAS DATA#############################
## The script assumes that data are in the following format:
## GWAS summary data is in gcta format (SNP, A1, A2, frq, b, se, p, N)
## eQTL summary data is in SMR format (SNP, Chr, BP, A1, A2, Freq, Probe, Probe_Chr, Probe_bp, Gene, Orientation, b, SE, p)
##########################################################################################

# Read in the list of enselmbl ids for genes to be analysed
genelist=read.table("./ENSG_IDs_for_genes_near_ACE.txt", header=T, stringsAsFactor=F, sep="\t")
gene_symbol=genelist$Gene.name
gene_ENSG_ID=genelist$Gene.stable.ID

# Parameters for schizophrenia data
path1="./"
file1="gwas_summary_data"
type1="cc" ## Type of trait - cc for case/ctrl
s1= 0.3862 # proportion of cases. For SCZ data this is 0.3862 (40,675 cases, 64,643 controls)
gwas_data=read.table(paste(path1,file1,sep=""), header=T, stringsAsFactor=F)

# Results will be save in this object
coloc_exp_results=c()

# Perform coloc analysis between eqtl and scz data for each gene
for (x in 1:nrow(genelist)) {
    g=gene_symbol[x]
    path2="./smr_query_files/"
    file2=paste(gene_ENSG_ID[x],"_",gene_symbol[x],"_eqtlgen_summary.txt", sep="")
    type2="quant" # quantitative trait

    # Read in eqtl
    d2=read.table(paste(path2,file2,sep=""), header=T, stringsAsFactor=F)
    d1=gwas_data
    
    ## Extract SNPs present in both datasets
    snps = intersect(d1$SNP,d2$SNP)

    ## Harmonise SNP order in the 2 datasets
    rownames(d1)=d1$SNP
    rownames(d2)=d2$SNP
    d1=d1[snps,]
    d2=d2[snps,]

    # assign paramter values and create list of input parameters for dataset 1
    MAF1=ifelse(d1$frq<0.5,d1$frq, 1-d1$frq)
    pvalues1=d1$p
    beta1=d1$b
    varbeta1=d1$se^2
    N1=d1$N
    list1=list("snps"=snps,"MAF"=MAF1,"beta"=beta1,"varbeta"=varbeta1,"pvalues"=pvalues1,"N"=N1, "type"=type1, "s"=s1)

    # assign paramter values and create list of input parameters for dataset 2
    MAF2=ifelse(d2$Freq<0.5,d2$Freq, 1-d2$Freq)
    pvalues2=d2$p
    beta2=d2$b
    varbeta2=d2$SE^2
    N2=31684
    ## sdY estimation done by coloc function if unknown
    list2=list("snps"=snps,"MAF"=MAF2,"beta"=beta2,"varbeta"=varbeta2,"pvalues"=pvalues2,"N"=N2, "type"=type2)

    res=coloc.abf(list1,list2, p1 = 1e-04, p2 = 1e-04,p12 = 1e-05)

    # save coloc output to results object
    coloc_exp_results$gene[x]=gene_symbol[x]
    coloc_exp_results$PP.H4.abf[x]=res$summary["PP.H4.abf"]

}

coloc_exp_results=as.data.frame(coloc_exp_scz_results)
##########################################################################################
