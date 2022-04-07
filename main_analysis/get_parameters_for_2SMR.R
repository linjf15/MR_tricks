library(TwoSampleMR)

# Calculating instrument F-statistic for Two-sample MR using beta and se
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5446088/
# F = b^2/se^2
get_fval_bowden_method <- function(b, se)
{
    Fval <- b^2/se^2
    return(Fval)
}

## Calculating F-statistic from QTL p-value and N
## This should give almost identical values to the Bowden method
get_fval_2SMR <- function(p,n)
{
    Fval <- suppressWarnings(qf(p, 1, n - 1, low = FALSE))
    return(Fval)
}

## Use function get_r_from_pn(p=d$p.value, n=d$N) to get r
get_rsq_2SMR <- function(p,n) {
    r <- get_r_from_pn(p, n)
    rsq <- r^2
    return(rsq)
}

## Calculate power
## The power is estimated using the method by Burgess et al.
## https://www.ncbi.nlm.nih.gov/pubmed/24608958
## and can also be estimated using the online power calculator:
## https://sb452.shinyapps.io/power/

# ratio is 1:X (case:ctrl)
# n=sample size,
# rsq is proportion of variance in exposure explained by SNP

calc_power <- function(b1, n, ratio, rsq) {
    power <- pnorm(sqrt(n*rsq*(ratio/(1+ratio))*(1/(1+ratio)))*b1-qnorm(1-sig/2))
    return(power)
}
