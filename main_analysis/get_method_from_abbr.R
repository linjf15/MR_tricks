get_method_from_abbr <- function(mr_method = c("Wald ratio","Inverse variance weighted"))
{
  requireNamespace("stringr",quietly = TRUE)
  mr_method_list <- lapply(1:length(mr_method), function(i)
    {
    mr_method_abbr_i <- trimws(tolower(mr_method[i]))
    if(grepl("wr|wald ratio", 
             mr_method_abbr_i)) return("Wald ratio")
    if(grepl("ivw|inverse variance weighted",
             mr_method_abbr_i)) return("Inverse variance weighted")
    if(grepl("mre|egger|mr egger", 
             mr_method_abbr_i)) return("MR Egger")
    if(grepl("wm|weighted median",
             mr_method_abbr_i)) return("Weighted median")
    if(grepl("raps|robust adjusted profile score",
             mr_method_abbr_i)) return("Robust adjusted profile score (RAPS)")
    if(grepl("presso",
             mr_method_abbr_i)) return("MR PRESSO")
  })
  unique(unlist(mr_method_list))
}
