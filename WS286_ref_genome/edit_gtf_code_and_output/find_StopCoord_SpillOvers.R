library(tidyverse)
library(magrittr)
library(readxl)
library(writexl)

ws286 <- read_xlsx("./ws286_new.xlsx", col_names = T)
# ws286 <- read_xlsx("../ext_3pUtr_code_and outputs/ws286_ext_Gene3pUtrTransExon_50bp.xlsx", col_names = T)

# find genes whose stop position spill over into other genes' regions
ws286$seqname %<>% as.factor() %>% fct_relevel("I", "II", "III", "IV", "V", "X", "MtDNA")
ws286 %<>% arrange(seqname, gene_id, start)
ws286$StopCoordSpillOver <- NA
ws286$SpillIntoWhichGenes <-NA

ws286_find_spillovers_I <- ws286 %>% filter(feature == "gene" & seqname == "I")
ws286_find_spillovers_II <- ws286 %>% filter(feature == "gene" & seqname == "II")
ws286_find_spillovers_III <- ws286 %>% filter(feature == "gene" & seqname == "III")
ws286_find_spillovers_IV <- ws286 %>% filter(feature == "gene" & seqname == "IV")
ws286_find_spillovers_V <- ws286 %>% filter(feature == "gene" & seqname == "V")
ws286_find_spillovers_X <- ws286 %>% filter(feature == "gene" & seqname == "X")
ws286_find_spillovers_MtDNA <- ws286 %>% filter(feature == "gene" & seqname == "MtDNA")

find_spillovers <- function(df) {
  for (i in 1:nrow(df)) {
    if(df$strand[i] == "+") {
      if_spill <- apply(df[-i, c("start", "end")], MARGIN = 1, function(vec) between(df$end[i], vec[1], vec[2]))
      df$StopCoordSpillOver[i] <- any(if_spill)
      df$SpillIntoWhichGenes[i] <- df[-i, "gene_id"] %>% .[which(if_spill), ] %>% unlist() %>% paste(collapse = ",")
    } else {
      if_spill <- apply(df[-i, c("start", "end")], MARGIN = 1, function(vec) between(df$start[i], vec[1], vec[2]))
      df$StopCoordSpillOver[i] <- any(if_spill)
      df$SpillIntoWhichGenes[i] <- df[-i, "gene_id"] %>% .[which(if_spill), ] %>% unlist() %>% paste(collapse = ",")
    }
  }
 return(df) 
}

ws286_find_spillovers_I %<>% find_spillovers()
ws286_find_spillovers_II %<>% find_spillovers()
ws286_find_spillovers_III %<>% find_spillovers()
ws286_find_spillovers_IV %<>% find_spillovers()
ws286_find_spillovers_V %<>% find_spillovers()
ws286_find_spillovers_X %<>% find_spillovers()
ws286_find_spillovers_MtDNA %<>% find_spillovers()

ws286_find_spillovers <- rbind(ws286_find_spillovers_I,
                               ws286_find_spillovers_II,
                               ws286_find_spillovers_III,
                               ws286_find_spillovers_IV,
                               ws286_find_spillovers_V,
                               ws286_find_spillovers_X,
                               ws286_find_spillovers_MtDNA)

write_xlsx(ws286_find_spillovers, "ws286_find_spillovers.xlsx", col_names = T)
# write_xlsx(ws286_find_spillovers, "ws286_50bp_find_spillovers.xlsx", col_names = T)

# exclude cases where two genes are on different strands and one gene stops at the other gene's beginning
ws286_spillover <- read_xlsx("ws286_find_spillovers.xlsx", col_names = T)
ws286_spillover$more_than_one_spillover <- str_detect(ws286_spillover$SpillIntoWhichGenes, ",")
ws286_spillover$should_exclude <- NA

for(i in which(ws286_spillover$more_than_one_spillover == FALSE)) {
  spillover_gene <- ws286_spillover %>% filter(gene_id == ws286_spillover$SpillIntoWhichGenes[i])
  if_diff_strand <- (ws286_spillover$strand[i] != spillover_gene$strand)
  if (if_diff_strand) {
    if (ws286_spillover$strand[i] == "+") {
      if_end_equals_start <- (ws286_spillover$end[i] == spillover_gene$end)
    } else {
      if_end_equals_start <- (ws286_spillover$start[i] == spillover_gene$start)
    }
    ws286_spillover$should_exclude[i] <- if_end_equals_start
  } else {
    ws286_spillover$should_exclude[i] <- FALSE
  }
}

ws286_spillover$should_exclude %>% which() %>% length() # 725

for(i in which(ws286_spillover$more_than_one_spillover == TRUE)) {
  spillover_gene <- ws286_spillover %>% filter(gene_id %in% unlist(str_split(ws286_spillover$SpillIntoWhichGenes[i], ",")))
  if_diff_strand <- (ws286_spillover$strand[i] != spillover_gene$strand) %>% all()
  if (if_diff_strand) {
    if (ws286_spillover$strand[i] == "+") {
      if_end_equals_start <- (ws286_spillover$end[i] == spillover_gene$end) %>% all()
    } else {
      if_end_equals_start <- (ws286_spillover$start[i] == spillover_gene$start) %>% all()
    }
    ws286_spillover$should_exclude[i] <- if_end_equals_start
  } else {
    ws286_spillover$should_exclude[i] <- FALSE
  }
}

ws286_spillover$should_exclude %>% which() %>% length() # still 725

ws286_spillover$should_exclude[is.na(ws286_spillover$should_exclude)] <- "FALSE"
ws286_spillover$should_exclude %<>% as.logical()
ws286_spillover$should_exclude %>% which() %>% length() # still 725

write_xlsx(ws286_spillover, "ws286_find_spillovers_ShouldYouExclSome.xlsx", col_names = T)

