get_name_from_abbr <- function(mr_method = c("Wald ratio","Inverse variance weighted"))
{
  requireNamespace("stringr",quietly = TRUE)
  mr_method_list <- c()
  for (i in 1:length(mr_method)) {
    if(str_detect(mr_method[i], "wr|WR|Wald ratio|wald ratio"))
    {
      mr_method_list <- c(mr_method_list,"Wald ratio")
    }
    if(str_detect(mr_method[i], "ivw|IVW|Inverse variance weighted|inverse variance weighted"))
    {
      mr_method_list <- c(mr_method_list,"Inverse variance weighted")
    }
    if(str_detect(mr_method[i], "mre|MRE|Egger|egger"))
    {
      mr_method_list <- c(mr_method_list,"MR Egger")
    }
    if(str_detect(mr_method[i], "wm|WM|Weighted median|weighted median"))
    {
      mr_method_list <- c(mr_method_list,"Weighted median")
    }
    if(str_detect(mr_method[i], "raps|RAPS|Robust adjusted profile score|robust adjusted profile score"))
    {
      mr_method_list <- c(mr_method_list,"Robust adjusted profile score (RAPS)")
    }
    if(str_detect(mr_method[i], "presso|PRESSO"))
    {
      mr_method_list <- c(mr_method_list,"MR PRESSO")
    }
  }
  mr_method_list
}
