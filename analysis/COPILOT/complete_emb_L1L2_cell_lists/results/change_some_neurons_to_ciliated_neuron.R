library(tidyverse)
library(magrittr)
library(readxl)
library(writexl)

cil <- read_tsv("worm_ciliated_neurons.txt", col_names = F)
pat <- paste0(cil$X1, collapse = "|")

emb <- read_xlsx("at_hatch_558_cells.xlsx")
temp <- emb[str_detect(emb$cell_common_name, pat), ]
temp <- temp[str_detect(temp$cell_common_name, "so|sh", negate = T), ]
temp$cell_type <- "ciliated_neuron"
emb %<>% filter(!(cell_common_name %in% temp$cell_common_name))
emb %<>% rbind(temp)
write_xlsx(emb, "at_hatch_558_cells.xlsx")

l1 <- read_xlsx("6h_fed_l1.xlsx")
temp <- l1[str_detect(l1$cell_common_name, pat), ]
temp <- temp[str_detect(temp$cell_common_name, "so|sh", negate = T), ]
temp$cell_type <- "ciliated_neuron"
l1 %<>% filter(!(cell_common_name %in% temp$cell_common_name))
l1 %<>% rbind(temp)
write_xlsx(l1, "6h_fed_l1.xlsx")
