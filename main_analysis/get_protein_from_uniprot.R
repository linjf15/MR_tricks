get_protein_from_uniprot <- function(UniProtID,local = TRUE,ref_panel = "UniProtRefPanel.xlsx",
                                     col_uniprot = "Entry",col_protein = "Protein")
{
  if(!local)
  {
    requireNamespace("UniprotR", quietly = TRUE)
    uniprot_query <- try(GetProteinAnnontate(trimws(UniProtID),
                                             columns = "protein names"))
    if("try-error" %in% class(uniprot_query)) 
    {
      print(paste0("Researching for ", UniProtID))
      uniprot_query <- try(GetProteinAnnontate(trimws(UniProtID),
                                               columns = "protein names"))
      if("try-error" %in% class(uniprot_query)) 
      {
        print(paste0("Researching for ", UniProtID, "Failed"))
        return(paste0("https://www.uniprot.org/uniprot/",toupper(trimws(UniProtID))))
      }
      else
      {
        return(uniprot_query)
      }
    }
    else
    {
      return(uniprot_query)
    }
  }
  
  else
  {
    ref_panel <- readxl::read_excel(ref_panel)
    ref_panel_query <- ref_panel[ref_panel[[col_uniprot]] == UniProtID,][[col_protein]]
    if(length(ref_panel_query)==1) return(ref_panel_query)
    else return(paste0("https://www.uniprot.org/uniprot/",toupper(trimws(UniProtID))))
  }
}
