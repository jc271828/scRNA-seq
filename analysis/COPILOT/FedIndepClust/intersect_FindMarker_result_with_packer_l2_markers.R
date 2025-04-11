library(tidyverse)
library(magrittr)
library(readxl)
library(writexl)

multiplesheets <- function(fname) {
  # getting info about all excel sheets
  sheets <- excel_sheets(fname)
  tibble <- lapply(sheets, function(x) read_xlsx(fname, sheet = x))
  data_frame <- lapply(tibble, as.data.frame)
  # assigning names to data frames
  names(data_frame) <- sheets
  # print data frame
  print(data_frame)
}

marker0.2 <- multiplesheets("~/../OneDrive - Duke University/lab/experiments/20220719_scRNASeq/analysis/COPILOT/FilRatio1TopPercVariesMinUmi1025MtThresh20/FedIndepClust/FedRes0.2MarkerGenes_20230525.xlsx")

# packer_emb_df <- read_xlsx("~/../OneDrive - Duke University/lab/COMMON_FILES/Packer2019_scRNAseq_embryo/Packer2019_emb_cell_type_signature_genes.xlsx")
# packer_l2_df <- read_xlsx("~/../OneDrive - Duke University/lab/COMMON_FILES/Packer2019_scRNAseq_embryo/Packer2019_L2_cell_type_signature_genes.xlsx")
# 
# packer_emb_df$exclude[is.na(packer_emb_df$exclude)] <- "NA"
# packer_l2_df$exclude[is.na(packer_l2_df$exclude)] <- "NA"
# 
# packer_emb_df$marker_genes %<>% paste(",", sep = "")
# packer_emb_df$exclude[packer_emb_df$exclude != "NA"] %<>% paste(",", sep = "")
# packer_l2_df$marker_genes %<>% paste(",", sep = "")
# packer_l2_df$exclude[packer_l2_df$exclude != "NA"] %<>% paste(",", sep = "")

# annot_cluster <- function(df) {
#   df %<>% select(-c("packer", "packer_not", "cao", "cao_not"))
#   
#   df$packer_emb <- NA
#   df$packer_emb_not <- NA
#   df$packer_l2 <- NA
#   df$packer_l2_not <- NA
#   
#   for (i in which(df$avg_log2FC > 0 & df$p_val_adj < 0.05)) {
#     df[i, "packer_emb"] <- packer_emb_df$cell_type[str_detect(packer_emb_df$marker_genes, paste(df$gene_name[i], ",", sep = ""))] %>% paste(collapse = "; ")
#     df[i, "packer_l2"] <- packer_l2_df$cell_type[str_detect(packer_l2_df$marker_genes, paste(df$gene_name[i], ",", sep = ""))] %>% paste(collapse = "; ")
#     
#     df[i, "packer_emb_not"] <- packer_emb_df$cell_type[str_detect(packer_emb_df$exclude, paste(df$gene_name[i], ",", sep = ""))] %>% paste(collapse = "; ")
#     df[i, "packer_l2_not"] <- packer_l2_df$cell_type[str_detect(packer_l2_df$exclude, paste(df$gene_name[i], ",", sep = ""))] %>% paste(collapse = "; ")
#   }
#   return(df)
# }

garnett_cao_l2_df <- read_xlsx("~/../OneDrive - Duke University/lab/COMMON_FILES/Cao2017_scRNAseq_early_L2s/GarnettTrainedMarkers/Cao2017_L2_cell_type_signature_genes_GarnettTrained.xlsx")
garnett_cao_l2_df$exclude[is.na(garnett_cao_l2_df$exclude)] <- "NA"
garnett_cao_l2_df$marker_genes %<>% paste(",", sep = "")
garnett_cao_l2_df$exclude[garnett_cao_l2_df$exclude != "NA"] %<>% paste(",", sep = "")

annot_cluster <- function(df) {
  df$garnett_cao_l2 <- NA

  for (i in which(df$avg_log2FC > 0 & df$p_val_adj < 0.05)) {
    df[i, "garnett_cao_l2"] <- garnett_cao_l2_df$cell_type[str_detect(garnett_cao_l2_df$marker_genes, paste(df$gene_name[i], ",", sep = ""))] %>% paste(collapse = "; ")
  }
  return(df)
}
  
for (j in 1:length(marker0.2)) {
  df <- marker0.2[[j]] %>% annot_cluster()
  marker0.2[[j]] <- df
}

temp <- marker0.2[[1]]

write_xlsx(marker0.2, "FedRes0.2MarkerGenes_20230531.xlsx")
