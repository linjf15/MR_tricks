```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

This is a modified copy from [Fang Li](https://hanruizhang.github.io/GWAS-eQTL-Colocalization/).

### 1. The purpose of GWAS-eQTL intergration
* Is the my variant an eQTL?
* Is the leading variant of the GWAS and eQTL signal the same?
* Is my GWAS association of interest driven by an eQTL that may indiciate a functinal mechanism?  

GWAS locus that colocalized with eQTL is one of the primary and scalable signal for functional follow-up analyses. 

### 2. Install R/RStudio and packages
* Install the latest version of [R](https://cran.r-project.org/) or [RStudio](https://www.rstudio.com/products/rstudio/download/).
* Install R pakage [locuscomparer](https://github.com/boxiangliu/locuscomparer). 
* Install R package ***coloc***:   
  + `if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager") `   
  + `BiocManager::install("snpStats")`
  + `install("coloc")`

### 3. Colocalization analysis using ***coloc***  
* Read sample data into R:  
      You can download the examples files:     [GWAS](https://github.com/fangli359/fangli359.github.io/blob/master/CAD_GWAS.txt) and [eQTL](https://github.com/fangli359/fangli359.github.io/blob/master/Artery_Coronary_v7_eQTL_PHACTR1.txt) datasets.  
  + `eqtl <- read.table(file="[path to]/Artery_Coronary_v7_eQTL_PHACTR1.txt", header=T, as.is=T)`   
  + `gwas <- read.table(file="[path to]/CAD_GWAS.txt", header=T, as.is=T)`
   
* Merge gwas and eqtl data sets by only shared "rs _ id":  
  + `input <- merge(eqtl, gwas, by="rs_id", all=FALSE, suffixes=c("_eqtl","gwas")`  

    Optinal: provide suffix to differentiate data source from gwas or eqtl.   

* Run ***coloc*** using [*coloc.abf()*](https://cran.r-project.org/web/packages/coloc/vignettes/vignette.html) fuction:  
  + `dataset1 <- list(pvalues=input$pval_nominal_gwas, type="cc", S=0.33, N=nrow(gwas))`
  + `dataset2 <- list(pvalues=input$pval_nominal_eqtl, type="quant", N=nrow(eqtl))`
  + `result <- coloc.abf(dataset1, dataset2, MAF=input$maf)`   

    Comments: [*coloc.abf()*](https://cran.r-project.org/web/packages/coloc/vignettes/vignette.html) function needs two named lists (gwas and eqtl) that contain p-values, the type of study ("cc" for case-control studies, "quant" for quantitative traits) and sample size(N). s= the proportion of samples are cases, when type="cc". It also needs the minor allele frequency. 
 
* Read out posterior probabilities for colocalization:  

       H0: neither trait has a genetic association in the region  
   
       H1/H2: only trait 1/trait 2 has a genetic association in the region  
   
       H3: both traits are associated, but with different causal variants  
   
       H4: both traits are associated and share a single causal variant  
   
     A posterior probability of ≥75%/ 80% is considered strong evidence of the eQTL-GWAS pair influencing both the expression and GWAS trait at a particular region.   
 
### 4. Visualization using *locuscomparer* 
* Define file names of the GWAS and eQTL data sets:  
  + `gwas_fn <- "CAD_GWAS.txt"`  
  + `eqtl_fn <- "Artery_Coronary_v7_eQTL_PHACTR1.txt"`  
  + `marker_col <- "rs_id"`  
  + `pval_col <- "pval_nominal"`    
  
* Run *locuscompare* to visualize:
  + `locuscompare(in_fn1=gwas_fn, in_fn2=eqtl_fn, title1="GWAS", title2="eQTL", marker_col1= marker_col, pval_col1=pval_col, marker_col2=marker_col, pval_col2=pval_col))`  
  
* Results output:   
   
![](https://raw.githubusercontent.com/fangli359/fangli359.github.io/master/rs9349379%20locus.png)
