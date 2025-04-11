library(tidyverse)
library(magrittr)
library(readxl)
library(writexl)

setwd(".")

#########################################################################################################

# cell_names <- read_xlsx("complete_cell_list_wormatlas_sulston.xlsx")
# cell_names$`Lineage Name` %<>% str_remove_all(" ") %>% str_remove("\\.")
# 
# cells0 <- cell_names$`Lineage Name` %>% unique()
# cells1 <- cells0 %>% sapply(FUN = function(x) substr(x, start = 1, stop = nchar(x) - 1)) %>% unique()
# cells2 <- cells1 %>% sapply(FUN = function(x) substr(x, start = 1, stop = nchar(x) - 1)) %>% unique()
# cells3 <- cells2 %>% sapply(FUN = function(x) substr(x, start = 1, stop = nchar(x) - 1)) %>% unique()
# cells4 <- cells3 %>% sapply(FUN = function(x) substr(x, start = 1, stop = nchar(x) - 1)) %>% unique()
# cells5 <- cells4 %>% sapply(FUN = function(x) substr(x, start = 1, stop = nchar(x) - 1)) %>% unique()
# cells6 <- cells5 %>% sapply(FUN = function(x) substr(x, start = 1, stop = nchar(x) - 1)) %>% unique()
# cells7 <- cells6 %>% sapply(FUN = function(x) substr(x, start = 1, stop = nchar(x) - 1)) %>% unique()
# cells8 <- cells7 %>% sapply(FUN = function(x) substr(x, start = 1, stop = nchar(x) - 1)) %>% unique()
# cells9 <- cells8 %>% sapply(FUN = function(x) substr(x, start = 1, stop = nchar(x) - 1)) %>% unique()
# cells10 <- cells9 %>% sapply(FUN = function(x) substr(x, start = 1, stop = nchar(x) - 1)) %>% unique()
# cells11 <- cells10 %>% sapply(FUN = function(x) substr(x, start = 1, stop = nchar(x) - 1)) %>% unique()
# cells12 <- cells11 %>% sapply(FUN = function(x) substr(x, start = 1, stop = nchar(x) - 1)) %>% unique()
# cells13 <- cells12 %>% sapply(FUN = function(x) substr(x, start = 1, stop = nchar(x) - 1)) %>% unique()
# cells14 <- cells13 %>% sapply(FUN = function(x) substr(x, start = 1, stop = nchar(x) - 1)) %>% unique()
# cells15 <- cells14 %>% sapply(FUN = function(x) substr(x, start = 1, stop = nchar(x) - 1)) %>% unique()
# cells16 <- cells15 %>% sapply(FUN = function(x) substr(x, start = 1, stop = nchar(x) - 1)) %>% unique()
# cells17 <- cells16 %>% sapply(FUN = function(x) substr(x, start = 1, stop = nchar(x) - 1)) %>% unique()
# 
# all_possible_cells <-
#   c(cells0, cells1, cells2, cells3, cells4, cells5, cells6, cells7, cells8,
#     cells9, cells10, cells11, cells12, cells13, cells14, cells15, cells16, cells17) %>%
#   unique() %>%
#   sort()
# 
# all_possible_cells <- all_possible_cells[-1]
# 
# all_possible_cells %>% as.data.frame() %>% write_tsv("parameters_for_octoparse.txt", col_names = F)

#########################################################################################################

crawl_wormweb <- read_xlsx("wormweb_crawl_result.xlsx")
crawl_wormweb <- crawl_wormweb[!(duplicated(crawl_wormweb)), ]
crawl_wormweb %<>% arrange(Linealname, Field)

colnames(crawl_wormweb) <- c("field", "value", "cell")

crawl_wormweb_wide <- crawl_wormweb %>% pivot_wider(names_from = field, values_from = value, id_cols = cell)

crawl_wormweb_wide <- crawl_wormweb_wide[, -6]
colnames(crawl_wormweb_wide) <- c("cell_lineage_name", "birth_time", "cell_type", "cell_common_name", "death_time")

crawl_wormweb_wide$birth_time %<>% str_remove_all(" min")
crawl_wormweb_wide$death_time %<>% str_remove_all(" min")

crawl_wormweb_wide$birth_time %<>% as.numeric()
crawl_wormweb_wide$death_time %<>% as.numeric()

colnames(crawl_wormweb_wide)[2] <- "birth_time_min"
colnames(crawl_wormweb_wide)[5] <- "death_time_min"

crawl_wormweb_wide %<>% arrange(birth_time_min, cell_lineage_name)

crawl_wormweb_wide$born_in <- NA
crawl_wormweb_wide$born_in[crawl_wormweb_wide$birth_time_min < 800] <- "embryo"
crawl_wormweb_wide$born_in[crawl_wormweb_wide$birth_time_min > 800 & crawl_wormweb_wide$birth_time_min < 1120] <- "early_L1"
crawl_wormweb_wide$born_in[crawl_wormweb_wide$birth_time_min > 1120 & crawl_wormweb_wide$birth_time_min < 1440] <- "mid_L1"
crawl_wormweb_wide$born_in[crawl_wormweb_wide$birth_time_min > 1440 & crawl_wormweb_wide$birth_time_min < 1760] <- "late_L1"

crawl_wormweb_wide$born_in[crawl_wormweb_wide$birth_time_min > 1760 & crawl_wormweb_wide$birth_time_min < 1940] <- "early_L2"
crawl_wormweb_wide$born_in[crawl_wormweb_wide$birth_time_min > 1940 & crawl_wormweb_wide$birth_time_min < 2120] <- "mid_L2"
crawl_wormweb_wide$born_in[crawl_wormweb_wide$birth_time_min > 2120 & crawl_wormweb_wide$birth_time_min < 2300] <- "late_L2"

crawl_wormweb_wide$born_in[crawl_wormweb_wide$birth_time_min > 2300 & crawl_wormweb_wide$birth_time_min < 2480] <- "early_L3"
crawl_wormweb_wide$born_in[crawl_wormweb_wide$birth_time_min > 2480 & crawl_wormweb_wide$birth_time_min < 2660] <- "mid_L3"
crawl_wormweb_wide$born_in[crawl_wormweb_wide$birth_time_min > 2660 & crawl_wormweb_wide$birth_time_min < 2840] <- "late_L3"

crawl_wormweb_wide$born_in[crawl_wormweb_wide$birth_time_min > 2840 & crawl_wormweb_wide$birth_time_min < 3060] <- "early_L4"
crawl_wormweb_wide$born_in[crawl_wormweb_wide$birth_time_min > 3060 & crawl_wormweb_wide$birth_time_min < 3280] <- "mid_L4"
crawl_wormweb_wide$born_in[crawl_wormweb_wide$birth_time_min > 3280 & crawl_wormweb_wide$birth_time_min < 3500] <- "late_L4"

crawl_wormweb_wide$born_in[crawl_wormweb_wide$birth_time_min > 3500 ] <- "adult"

#########################################################################################################

crawl_wormweb_wide$note <- NA
on_the_edge <- crawl_wormweb_wide[is.na(crawl_wormweb_wide$born_in), ]

on_the_edge[1:24, "born_in"] <- "early_L2"
on_the_edge[1:24, "note"] <- "HypFIG1.jpg WormAtlas"

on_the_edge[27:28, "born_in"] <- "late_L2"
on_the_edge[27:28, "note"] <- "HypFIG1.jpg WormAtlas"

on_the_edge[29:68, "born_in"] <- "early_L3"
on_the_edge[29:68, "note"] <- "HypFIG1.jpg WormAtlas"

on_the_edge[79:90, "born_in"] <- "late_L3"
on_the_edge[79:90, "note"] <- "HypFIG1.jpg WormAtlas"

on_the_edge[119:122, "born_in"] <- "early_L4"
on_the_edge[119:122, "note"] <- "HypFIG1.jpg WormAtlas"

on_the_edge[139:174, "born_in"] <- "early_L4"
on_the_edge[139:174, "note"] <- "HypFIG1.jpg WormAtlas"

on_the_edge[25:26, "born_in"] <- "early_L2"
on_the_edge[25:26, "note"] <- "QGKlineages.jpg WormAtlas"

on_the_edge[69:74, "born_in"] <- "early_L3"
on_the_edge[69:74, "note"] <- "Z1Z4lineages.jpg WormAtlas"

on_the_edge[75:78, "born_in"] <- "mid_L3"
on_the_edge[75:78, "note"] <- "Z1Z4lineages.jpg WormAtlas"

on_the_edge[91:94, "born_in"] <- "mid_L3"
on_the_edge[91:94, "note"] <- "Z1Z4lineages.jpg WormAtlas"

on_the_edge[95:102, "born_in"] <- "late_L3"
on_the_edge[95:102, "note"] <- "Z1Z4lineages.jpg WormAtlas"

on_the_edge[103:114, "born_in"] <- "mid_L3"
on_the_edge[103:114, "note"] <- "Z1Z4lineages.jpg WormAtlas"

on_the_edge[115:118, "born_in"] <- "late_L3"
on_the_edge[115:118, "note"] <- "Z1Z4lineages.jpg WormAtlas"

on_the_edge[175:182, "born_in"] <- "late_L3"
on_the_edge[175:182, "note"] <- "Z1Z4lineages.jpg WormAtlas"

on_the_edge[183:210, "born_in"] <- "late_L3"
on_the_edge[183:210, "note"] <- "Z1Z4lineages.jpg WormAtlas"

on_the_edge[123:138, "born_in"] <- "late_L3"
on_the_edge[123:138, "note"] <- "NonstriatedMuscle_VulvalMuscles WormAtlas"

on_the_edge[211, "born_in"] <- "embryo"
on_the_edge[211, "birth_time_min"] <- 0

#########################################################################################################

crawl_wormweb_wide <-
  rbind(crawl_wormweb_wide[!(is.na(crawl_wormweb_wide$born_in)), ],
        on_the_edge) %>%
  arrange(birth_time_min, cell_lineage_name)

crawl_wormweb_wide <- crawl_wormweb_wide[is.na(crawl_wormweb_wide$death_time_min), ]

crawl_wormweb_wide %<>% select(-c("death_time_min"))

#########################################################################################################
crawl_wormweb_wide_raw <- crawl_wormweb_wide

crawl_wormweb_wide %<>% arrange(birth_time_min, cell_lineage_name)

another_emb_cell_list <- read_xlsx("JC_manually_curated_cell_list.xlsx") %>% filter(born_in == "embryo")

born_in_emb <-
  crawl_wormweb_wide[crawl_wormweb_wide$birth_time_min < 800, ] %>%
  rbind(another_emb_cell_list) %>%
  unique() %>%
  arrange(birth_time_min, cell_lineage_name)

born_in_emb$is_parent <- NA

for (i in 1:nrow(born_in_emb)) {
  temp1 <-
    born_in_emb$cell_lineage_name[-i] %>%
    str_detect(pattern = paste("^", born_in_emb$cell_lineage_name[i], sep = "")) %>%
    any()

  # temp2 <-
  #   born_in_emb$cell_lineage_name[-i] %>%
  #   str_detect(pattern = paste("^", born_in_emb$cell_common_name[i], sep = "")) %>%
  #   any()
  
  temp3 <-
    born_in_emb$cell_common_name[-i] %>%
    str_detect(pattern = paste("^", born_in_emb$cell_lineage_name[i], sep = "")) %>%
    any()
  # 
  # temp4 <-
  #   born_in_emb$cell_common_name[-i] %>%
  #   str_detect(pattern = paste("^", born_in_emb$cell_common_name[i], sep = "")) %>%
  #   any()
  
  born_in_emb$is_parent[i] <- any(temp1, temp3)
}

terminal <- born_in_emb[born_in_emb$is_parent == FALSE, ]

terminal %<>% arrange(birth_time_min, cell_lineage_name)
terminal <- terminal[-1, ]

#########################################################################################################
terminal %<>% arrange(cell_lineage_name)

terminal[str_detect(terminal$cell_common_name, "^int"), ] %>% arrange(cell_common_name) # correct, 20 total (refer to IntFIG3.jpg WormAtlas)

# refer to the print-out of embryonic lineage
terminal[str_detect(terminal$cell_lineage_name, "^ABalaa|^ABalap"), ] %>% arrange(cell_lineage_name) # correct, 42 total
terminal[str_detect(terminal$cell_lineage_name, "^ABalpa|^ABalpp"), ] %>% arrange(cell_lineage_name) # correct, 49 total
terminal[str_detect(terminal$cell_lineage_name, "^ABaraa|^ABarap"), ] %>% arrange(cell_lineage_name) # correct, 55 total
terminal[str_detect(terminal$cell_lineage_name, "^ABarpa|^ABarpp"), ] %>% arrange(cell_lineage_name) # correct, 37 total
terminal[str_detect(terminal$cell_lineage_name, "^ABplaa|^ABplap"), ] %>% arrange(cell_lineage_name) # correct, 43 total
terminal[str_detect(terminal$cell_lineage_name, "^ABplpa|^ABplpp"), ] %>% arrange(cell_lineage_name) # correct, 57 total
terminal[str_detect(terminal$cell_lineage_name, "^ABpraa|^ABprap"), ] %>% arrange(cell_lineage_name) # correct, 47 total
terminal[str_detect(terminal$cell_lineage_name, "^ABprpa|^ABprpp"), ] %>% arrange(cell_lineage_name) # correct, 59 total
terminal[str_detect(terminal$cell_lineage_name, "^MS"), ] %>% arrange(cell_lineage_name) # correct, 80 total
terminal[str_detect(terminal$cell_lineage_name, "^E"), ] %>% arrange(cell_lineage_name) # correct, 20 total
terminal[str_detect(terminal$cell_lineage_name, "^C"), ] %>% arrange(cell_lineage_name) # correct, 47 total

terminal[str_detect(terminal$cell_lineage_name, "^D"), ] %>% arrange(cell_lineage_name) # missing Dpapp, should be 20 total
terminal[str_detect(terminal$cell_lineage_name, "Z2"), ] %>% arrange(cell_lineage_name) # missing
terminal[str_detect(terminal$cell_lineage_name, "Z3"), ] %>% arrange(cell_lineage_name) # missing

#########################################################################################################

terminal <-
  terminal[, -7] %>%
  rbind(data.frame(cell_lineage_name = c("Dpapp", "P4p", "P4a"),
                   birth_time_min = c(277, 140, 140),
                   cell_type = c("muscle", "repro", "repro"),
                   cell_common_name = c("mu bod", "Z2", "Z3"),
                   born_in = rep("embryo", 3),
                   note = NA))

terminal %<>% arrange(birth_time_min, cell_lineage_name)

at_hatch <- terminal

write_xlsx(at_hatch, "at_hatch_558_cells.xlsx")

#########################################################################################################

at_hatch <- read_xlsx("at_hatch_558_cells.xlsx")

anL1L2_cell_list <- read_xlsx("JC_manually_curated_cell_list.xlsx") %>% filter(birth_time_min >= 800 & birth_time_min <= 2300) %>% arrange(cell_lineage_name)

born_in_L1L2 <-
  crawl_wormweb_wide[crawl_wormweb_wide$birth_time_min >= 800 & crawl_wormweb_wide$birth_time_min <= 2300, ] %>%
  rbind(anL1L2_cell_list) %>%
  unique() %>%
  arrange(cell_lineage_name)

write_xlsx(born_in_L1L2, "L1L2_AllCells_Toedit.xlsx")

# Then manually edit this spreadsheet

L1L2AllCells <- read_xlsx("L1L2_AllCells_Toedit.xlsx") %>% arrange(birth_time_min, cell_lineage_name)

write_xlsx(L1L2AllCells, "L1L2_AllCells_Jc_Edited.xlsx")
