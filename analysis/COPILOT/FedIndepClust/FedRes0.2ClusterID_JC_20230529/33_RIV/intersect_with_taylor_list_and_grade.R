library(tidyverse)
library(magrittr)
library(readxl)
library(writexl)

taylor <- read_xlsx("~/../OneDrive - Duke University/lab/COMMON_FILES/Taylor2021_neuron_scRNASeq_L4/taylor2019L4_marker_genes.xlsx")
taylor <- taylor[, 1:4]
taylor %<>% filter(type %in% c("neuron", "pha_neuron"))
taylor$marker_genes %<>% paste(",", sep = "")

taylor %<>% select(-c("exclude"))
# taylor$exclude[is.na(taylor$exclude)] <- "NA"
# taylor$exclude[taylor$exclude != "NA"] %<>% paste(",", sep = "")

cluster <- read_xlsx("~/../OneDrive/Downloads/FedRes0.2MarkerGenes_20230531.xlsx", sheet = "cluster33.markers")

cluster$taylor_l4 <- NA
for (i in which(cluster$avg_log2FC > 0 & cluster$p_val_adj < 0.05)) {
  cluster[i, "taylor_l4"] <- taylor$cell_type[str_detect(taylor$marker_genes, paste(cluster$gene_name[i], ",", sep = ""))] %>% paste(collapse = "; ")
}

cluster %<>% select(-c("packer_emb", "packer_emb_not", "packer_l2", "packer_l2_not", "garnett_cao_l2"))

write_xlsx(cluster, "cluster33_TaylorL4Markers.xlsx")

##############################################################################
# grading with Taylor list

taylor <- read_xlsx("~/../OneDrive - Duke University/lab/COMMON_FILES/Taylor2021_neuron_scRNASeq_L4/taylor2019L4_marker_genes.xlsx")
taylor <- taylor[, 1:4]
taylor %<>% filter(type %in% c("neuron", "pha_neuron"))
taylor %<>% select(-c("exclude"))
# important to have a following white space!
taylor$marker_genes %<>% paste(", ", sep = "")

# cluster <- read_xlsx("~/../OneDrive/Downloads/cluster28_0_TaylorL4Markers.xlsx")

cluster %<>% filter(p_val_adj < 0.05 & avg_log2FC > 0)

possible_cell_types <- cluster$taylor_l4 %>% str_split("; ") %>% unlist() %>% unique()
possible_cell_types <- possible_cell_types[which(possible_cell_types != "NA")]
possible_cell_types <- possible_cell_types[which(nchar(possible_cell_types) != 0)]
possible_cell_types

temp <- data.frame(cell_type = NA,
                   taylor_l4_marker = NA,
                   if_significant_by_seurat_l2 = NA,
                   prop_present = NA)

for (i in 1:length(possible_cell_types)) {
  cell_type_spe_marker <- 
    taylor %>%
    filter(cell_type == possible_cell_types[i]) %>%
    .$marker_genes %>%
    str_split(" ") %>%
    unlist()
  cell_type_spe_marker <- cell_type_spe_marker[nchar(cell_type_spe_marker) != 0]
  
  temp[(nrow(temp) + 1):(nrow(temp) + length(cell_type_spe_marker)), "cell_type"] <- possible_cell_types[i]
  temp[(nrow(temp) + 1 - length(cell_type_spe_marker)):nrow(temp), "taylor_l4_marker"] <- cell_type_spe_marker %>% str_remove(",")
  temp[(nrow(temp) + 1 - length(cell_type_spe_marker)):nrow(temp), "if_significant_by_seurat_l2"] <- sapply(X = cell_type_spe_marker,
                                                                                                            FUN = function(x) any(str_detect(paste(cluster$gene_name, ",", sep = ""), x)))
  temp[(nrow(temp) + 1 - length(cell_type_spe_marker)):nrow(temp), "prop_present"] <- round(sum(temp[(nrow(temp) + 1 - length(cell_type_spe_marker)):nrow(temp), "if_significant_by_seurat_l2"])/length(cell_type_spe_marker), 2)
}

temp %<>% arrange(desc(prop_present))
temp <- temp[-nrow(temp), ]

write_xlsx(temp, "cluster33_GradingTaylorL4.xlsx")
