library(tidyverse)
library(magrittr)
library(readxl)
library(writexl)
library(ggrepel)
library(ggpubr)

got <- read_xlsx("../qc_tables/qc_per_identity_20231006.xlsx", sheet = "qc_per_tissue")
got <- got[str_detect(got$tissue, "_starv"), ]
got$tissue %<>% str_remove_all("_starv")

starv_total_num_cells <- sum(got$num_barcode)
expected <- read_xlsx("~/../Downloads/at_hatch_558_cells_CellNumPerTissue_20231006.xlsx")

# didn't get "Gblast" "Wblast"
setdiff(expected$tissue, got$tissue)
setdiff(got$tissue, expected$tissue)

compr <-
  merge(expected, got, by = "tissue", all = F) %>%
  mutate(expected_num_in_expt = round(cell_num/558 * starv_total_num_cells))

colnames(compr)[2:3] <- c("num_cell_per_animal", "annotated_num_in_expt")

compr$annotated_to_expected_ratio <- compr$annotated_num_in_expt/compr$expected_num_in_expt
write_xlsx(compr, "starv_cell_num_correlation_20231006.xlsx")

options(ggrepel.max.overlaps = Inf)
ggplot(compr, aes(x = expected_num_in_expt, y = annotated_num_in_expt)) +
  stat_function(color = "red",
                linewidth = 1,
                fun = function(x) x) +
  geom_point(size = 2) +
  geom_smooth(method = "lm") +
  stat_cor(method = "pearson") +
  geom_text_repel(label = compr$tissue, size = 4) +
  geom_text(mapping = aes(x = 17000, y = 16000), label = "y = x", color = "red", size = 5) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  theme_classic() +
  theme(aspect.ratio = 1) +
  labs(x = "(Proportionally) expected number of cells",
       y = "Annotated number of cells",
       title = "starv scRNASeq obtained ~ expected cell number correlation")

# model doesn’t explain much of variation of the data but it is significant (better than not having a model)
cor.test(compr$expected_num_in_expt, compr$annotated_num_in_expt)


