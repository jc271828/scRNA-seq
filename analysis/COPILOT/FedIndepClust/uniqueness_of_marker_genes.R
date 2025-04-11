library(magrittr)
library(tidyverse)
library(readxl)
library(writexl)

packer_l2_df <- read_xlsx("~/../OneDrive - Duke University/lab/COMMON_FILES/Packer2019_scRNAseq_embryo/Packer2019_L2_cell_type_signature_genes.xlsx")
packer_l2_df$exclude[is.na(packer_l2_df$exclude)] <- "NA"
packer_l2_df$marker_genes %<>% paste(",", sep = "")
packer_l2_df$exclude[packer_l2_df$exclude != "NA"] %<>% paste(",", sep = "")

packer_emb_df <- read_xlsx("~/../OneDrive - Duke University/lab/COMMON_FILES/Packer2019_scRNAseq_embryo/Packer2019_Emb_cell_type_signature_genes.xlsx")
packer_emb_df$exclude[is.na(packer_emb_df$exclude)] <- "NA"
packer_emb_df$marker_genes %<>% paste(",", sep = "")
packer_emb_df$exclude[packer_emb_df$exclude != "NA"] %<>% paste(",", sep = "")

unique_to_which_cell_types <- function(gene_name) {
  gene <- paste(gene_name, sep = ",")
  packer_l2_df$cell_type[str_detect(packer_l2_df$marker_genes, gene)] %>%
    paste0(collapse = "; ") %>%
    paste(paste(gene_name, "identifies all these L2 cell types:"), ., sep = " ") %>%
    print()
  
  packer_emb_df$cell_type[str_detect(packer_emb_df$marker_genes, gene)] %>%
    paste0(collapse = "; ") %>%
    paste(paste(gene_name, "identifies all these Embryo cell types:"), ., sep = " ") %>%
    print()
}

"glr-2" %>% unique_to_which_cell_types()

 



