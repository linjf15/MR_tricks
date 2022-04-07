# Calculating instrument F-statistic using beta and se
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
