library(magrittr)
library(tidyverse)
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

# marker0.2 <- multiplesheets("~/../OneDrive - Duke University/lab/experiments/20220719_scRNASeq/analysis/COPILOT/FilRatio1TopPercVariesMinUmi1025MtThresh20/FedIndepClust/FedRes0.2MarkerGenes_20230531.xlsx")
marker0.2 <- multiplesheets("~/../OneDrive/Downloads/FedRes0.2MarkerGenesCluster20_20230605.xlsx")

####################################################################################################################################################################
# Grading Packer L2 marker identified cell types
packer_l2_df <- read_xlsx("~/../OneDrive - Duke University/lab/COMMON_FILES/Packer2019_scRNAseq_embryo/Packer2019_L2_cell_type_signature_genes.xlsx")
packer_l2_df$exclude[is.na(packer_l2_df$exclude)] <- "NA"

# important to have a following white space!
packer_l2_df$marker_genes %<>% paste(", ", sep = "")
packer_l2_df$exclude[packer_l2_df$exclude != "NA"] %<>% paste(", ", sep = "")

# for each sheet/cluster, get all name-present/identified/possible cell type names
# get all marker genes of that cell type from Packer's L2 marker list
# how many of each cell type's markers are identified by Seurat as significant (p.adj < 0.05 & log2FC > 1)
grade0.2 <- list()

for (j in 1:length(marker0.2)) {
  markers <- marker0.2[[j]] %>% filter(p_val_adj < 0.05 & avg_log2FC > 0)
  
  possible_cell_types <- markers$packer_l2 %>% str_split("; ") %>% unlist() %>% unique()
  possible_cell_types <- possible_cell_types[which(possible_cell_types != "NA")]
  possible_cell_types
  
  temp <- data.frame(cell_type = NA,
                     packer_l2_marker = NA,
                     if_significant_by_seurat_l2 = NA,
                     prop_present = NA)
  
  for (i in 1:length(possible_cell_types)) {
    cell_type_spe_marker <- 
      packer_l2_df %>%
      filter(cell_type == possible_cell_types[i]) %>%
      .$marker_genes %>%
      str_split(" ") %>%
      unlist()
    cell_type_spe_marker <- cell_type_spe_marker[nchar(cell_type_spe_marker) != 0]
    
    temp[(nrow(temp) + 1):(nrow(temp) + length(cell_type_spe_marker)), "cell_type"] <- possible_cell_types[i]
    temp[(nrow(temp) + 1 - length(cell_type_spe_marker)):nrow(temp), "packer_l2_marker"] <- cell_type_spe_marker %>% str_remove(",")
    temp[(nrow(temp) + 1 - length(cell_type_spe_marker)):nrow(temp), "if_significant_by_seurat_l2"] <- sapply(X = cell_type_spe_marker,
                                                                                                              FUN = function(x) any(str_detect(paste(markers$gene_name, ",", sep = ""), x)))
    temp[(nrow(temp) + 1 - length(cell_type_spe_marker)):nrow(temp), "prop_present"] <- round(sum(temp[(nrow(temp) + 1 - length(cell_type_spe_marker)):nrow(temp), "if_significant_by_seurat_l2"])/length(cell_type_spe_marker), 2)
  }
  
  impossible_cell_types <- markers %>% filter(packer_l2_not != "NA") %>% select(c("gene_name", "packer_l2_not"))
  impossible_cell_types$concat <- paste(impossible_cell_types$packer_l2_not, impossible_cell_types$gene_name, sep = " by ")
  
  temp %<>% arrange(desc(prop_present))
  temp <- temp[c(nrow(temp), 1:(nrow(temp) - 1)), ]
  
  temp[1, 1] <- paste("NOT these cell types:", paste0(impossible_cell_types$concat, collapse = " OR "), sep = " ")
  
  grade0.2[[j]] <- temp
}

names(grade0.2) <- names(marker0.2)

write_xlsx(grade0.2, "FedRes0.2MarkerGenesGradeL2Cluster20_20230605.xlsx")

####################################################################################################################################################################
# Do the same thing for Packer Embryo Markers
packer_emb_df <- read_xlsx("~/../OneDrive - Duke University/lab/COMMON_FILES/Packer2019_scRNAseq_embryo/Packer2019_Emb_cell_type_signature_genes.xlsx")
packer_emb_df$exclude[is.na(packer_emb_df$exclude)] <- "NA"

# important to have a following white space!
packer_emb_df$marker_genes %<>% paste(", ", sep = "")
packer_emb_df$exclude[packer_emb_df$exclude != "NA"] %<>% paste(", ", sep = "")

# for each sheet/cluster, get all name-present/identified/possible cell type names
# get all marker genes of that cell type from Packer's Emb marker list
# how many of each cell type's markers are identified by Seurat as significant (p.adj < 0.05 & log2FC > 1)
grade0.2 <- list()

for (j in 1:length(marker0.2)) {
  markers <- marker0.2[[j]] %>% filter(p_val_adj < 0.05 & avg_log2FC > 0)
  
  possible_cell_types <- markers$packer_emb %>% str_split("; ") %>% unlist() %>% unique()
  possible_cell_types <- possible_cell_types[which(possible_cell_types != "NA")]
  possible_cell_types
  
  temp <- data.frame(cell_type = NA,
                     packer_emb_marker = NA,
                     if_significant_by_seurat_emb = NA,
                     prop_present = NA)
  
  for (i in 1:length(possible_cell_types)) {
    cell_type_spe_marker <- 
      packer_emb_df %>%
      filter(cell_type == possible_cell_types[i]) %>%
      .$marker_genes %>%
      str_split(" ") %>%
      unlist()
    cell_type_spe_marker <- cell_type_spe_marker[nchar(cell_type_spe_marker) != 0]
    
    temp[(nrow(temp) + 1):(nrow(temp) + length(cell_type_spe_marker)), "cell_type"] <- possible_cell_types[i]
    temp[(nrow(temp) + 1 - length(cell_type_spe_marker)):nrow(temp), "packer_emb_marker"] <- cell_type_spe_marker %>% str_remove(",")
    temp[(nrow(temp) + 1 - length(cell_type_spe_marker)):nrow(temp), "if_significant_by_seurat_emb"] <- sapply(X = cell_type_spe_marker,
                                                                                                              FUN = function(x) any(str_detect(paste(markers$gene_name, ",", sep = ""), x)))
    temp[(nrow(temp) + 1 - length(cell_type_spe_marker)):nrow(temp), "prop_present"] <- round(sum(temp[(nrow(temp) + 1 - length(cell_type_spe_marker)):nrow(temp), "if_significant_by_seurat_emb"])/length(cell_type_spe_marker), 2)
  }
  
  impossible_cell_types <- markers %>% filter(packer_emb_not != "NA") %>% select(c("gene_name", "packer_emb_not"))
  impossible_cell_types$concat <- paste(impossible_cell_types$packer_emb_not, impossible_cell_types$gene_name, sep = " by ")
  
  temp %<>% arrange(desc(prop_present))
  temp <- temp[c(nrow(temp), 1:(nrow(temp) - 1)), ]
  
  temp[1, 1] <- paste("NOT these cell types:", paste0(impossible_cell_types$concat, collapse = " OR "), sep = " ")
  
  grade0.2[[j]] <- temp
}

names(grade0.2) <- names(marker0.2)

write_xlsx(grade0.2, "FedRes0.2MarkerGenesGradeEmbCluster20_20230605.xlsx")
####################################################################################################################################################################
# Do the same thing for Garnett trained Cao L2 Markers
garnett_cao_l2_df <- read_xlsx("~/../OneDrive - Duke University/lab/COMMON_FILES/Cao2017_scRNAseq_early_L2s/GarnettTrainedMarkers/Cao2017_L2_cell_type_signature_genes_GarnettTrained.xlsx")
garnett_cao_l2_df[c(4, 11, 27), "cell_type"] <- c("Intestinal_rectal muscle", "Am_Ph sheath cells", "flp-1_plus interneurons")
garnett_cao_l2_df$exclude[is.na(garnett_cao_l2_df$exclude)] <- "NA"

# important to have a following white space!
garnett_cao_l2_df$marker_genes %<>% paste(", ", sep = "")
garnett_cao_l2_df$exclude[garnett_cao_l2_df$exclude != "NA"] %<>% paste(", ", sep = "")

# for each sheet/cluster, get all name-present/identified/possible cell type names
# get all marker genes of that cell type from Packer's Emb marker list
# how many of each cell type's markers are identified by Seurat as significant (p.adj < 0.05 & log2FC > 1)
grade0.2 <- list()

for (j in 1:length(marker0.2)) {
  markers <- marker0.2[[j]] %>% filter(p_val_adj < 0.05 & avg_log2FC > 0)
  
  possible_cell_types <- markers$garnett_cao_l2 %>% str_split("; ") %>% unlist() %>% unique()
  possible_cell_types <- possible_cell_types[which(possible_cell_types != "NA")]
  possible_cell_types
  
  if (length(possible_cell_types) == 0) {
    temp <- data.frame(cell_type = "No cell type marker being significant",
                       garnett_cao_l2_marker = NA,
                       if_significant_by_seurat_emb = NA,
                       prop_present = NA)
  } else {
    temp <- data.frame(cell_type = NA,
                       garnett_cao_l2_marker = NA,
                       if_significant_by_seurat_emb = NA,
                       prop_present = NA)
    
    for (i in 1:length(possible_cell_types)) {
      cell_type_spe_marker <- 
        garnett_cao_l2_df %>%
        filter(cell_type == possible_cell_types[i]) %>%
        .$marker_genes %>%
        str_split(" ") %>%
        unlist()
      cell_type_spe_marker <- cell_type_spe_marker[nchar(cell_type_spe_marker) != 0]
      
      temp[(nrow(temp) + 1):(nrow(temp) + length(cell_type_spe_marker)), "cell_type"] <- possible_cell_types[i]
      temp[(nrow(temp) + 1 - length(cell_type_spe_marker)):nrow(temp), "garnett_cao_l2_marker"] <- cell_type_spe_marker %>% str_remove(",")
      temp[(nrow(temp) + 1 - length(cell_type_spe_marker)):nrow(temp), "if_significant_by_seurat_emb"] <- sapply(X = cell_type_spe_marker,
                                                                                                                 FUN = function(x) any(str_detect(paste(markers$gene_name, ",", sep = ""), x)))
      temp[(nrow(temp) + 1 - length(cell_type_spe_marker)):nrow(temp), "prop_present"] <- round(sum(temp[(nrow(temp) + 1 - length(cell_type_spe_marker)):nrow(temp), "if_significant_by_seurat_emb"])/length(cell_type_spe_marker), 2)
    }
    
    temp %<>% arrange(desc(prop_present))
    temp <- temp[1:(nrow(temp) - 1), ]
  }
  
  grade0.2[[j]] <- temp
}

names(grade0.2) <- names(marker0.2)

write_xlsx(grade0.2, "FedRes0.2MarkerGenesGradeGarnettL2Cluster20_20230605.xlsx")








