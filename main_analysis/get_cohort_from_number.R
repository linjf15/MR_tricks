get_cohort_from_number <- function(x){
  cohorts <- lapply(x, function(i)
    {
    switch(i, "1" = "IMSGC", "2" = "UK Biobank", "3" = "FinnGen", "meta" = "Pooled")
  })
  unlist(cohorts)
}
