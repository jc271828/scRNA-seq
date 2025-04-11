library(tidyverse)
library(magrittr)
library(readxl)
library(writexl)

emb <- read_xlsx("at_hatch_558_cells.xlsx")
l1 <- read_xlsx("6h_fed_l1.xlsx")

emb_num <- emb %>% group_by(cell_type) %>% reframe(cell_num = n()) %>% arrange(desc(cell_num), cell_type)
colnames(emb_num)[1] <- "tissue"

l1_num <- l1 %>% group_by(cell_type) %>% reframe(cell_num = n()) %>% arrange(desc(cell_num), cell_type)
colnames(l1_num)[1] <- "tissue"

write_xlsx(emb_num, "at_hatch_558_cells_CellNumPerTissue.xlsx")
write_xlsx(l1_num, "6h_fed_l1_CellNumPerTissue.xlsx")
