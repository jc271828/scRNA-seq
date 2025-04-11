library(tidyverse)
library(magrittr)
library(readxl)
library(writexl)

ws286 <- read_xlsx("./ws286_ext_Gene3pUtrTrans_250bp_20230422.xlsx", col_names = T)

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
  df[nchar(df$SpillIntoWhichGenes) == 0, "SpillIntoWhichGenes"] <- NA
  return(df) 
}

ws286_find_spillovers_I %<>% find_spillovers()
# write_xlsx(ws286_find_spillovers_I, "temp1.xlsx")
ws286_find_spillovers_II %<>% find_spillovers()
# write_xlsx(ws286_find_spillovers_II, "temp2.xlsx")
ws286_find_spillovers_III %<>% find_spillovers()
# write_xlsx(ws286_find_spillovers_III, "temp3.xlsx")
ws286_find_spillovers_IV %<>% find_spillovers()
# write_xlsx(ws286_find_spillovers_IV, "temp4.xlsx")
ws286_find_spillovers_V %<>% find_spillovers()
# write_xlsx(ws286_find_spillovers_V, "temp5.xlsx")
ws286_find_spillovers_X %<>% find_spillovers()
# write_xlsx(ws286_find_spillovers_X, "tempX.xlsx")
ws286_find_spillovers_MtDNA %<>% find_spillovers()
# write_xlsx(ws286_find_spillovers_MtDNA, "tempM.xlsx")

ws286_find_spillovers <- rbind(ws286_find_spillovers_I,
                               ws286_find_spillovers_II,
                               ws286_find_spillovers_III,
                               ws286_find_spillovers_IV,
                               ws286_find_spillovers_V,
                               ws286_find_spillovers_X,
                               ws286_find_spillovers_MtDNA)

# write_xlsx(ws286_find_spillovers, "ws286_find_spillovers_50bp.xlsx", col_names = T)

# exclude cases where two genes are on different strands and one gene stops at the other gene's beginning
ws286_find_spillovers$more_than_one_spillover <- str_detect(ws286_find_spillovers$SpillIntoWhichGenes, ",")
ws286_find_spillovers$should_exclude <- NA

for(i in which(ws286_find_spillovers$more_than_one_spillover == FALSE)) {
  spillover_gene <- ws286_find_spillovers %>% filter(gene_id == ws286_find_spillovers$SpillIntoWhichGenes[i])
  if_diff_strand <- (ws286_find_spillovers$strand[i] != spillover_gene$strand)
  if (if_diff_strand) {
    if (ws286_find_spillovers$strand[i] == "+") {
      if_end_equals_start <- (ws286_find_spillovers$end[i] == spillover_gene$end)
    } else {
      if_end_equals_start <- (ws286_find_spillovers$start[i] == spillover_gene$start)
    }
    ws286_find_spillovers$should_exclude[i] <- if_end_equals_start
  } else {
    ws286_find_spillovers$should_exclude[i] <- FALSE
  }
}

ws286_find_spillovers$should_exclude %>% which() %>% length() # what number is this?

for(i in which(ws286_find_spillovers$more_than_one_spillover == TRUE)) {
  spillover_gene <- ws286_find_spillovers %>% filter(gene_id %in% unlist(str_split(ws286_find_spillovers$SpillIntoWhichGenes[i], ",")))
  if_diff_strand <- (ws286_find_spillovers$strand[i] != spillover_gene$strand) %>% all()
  if (if_diff_strand) {
    if (ws286_find_spillovers$strand[i] == "+") {
      if_end_equals_start <- (ws286_find_spillovers$end[i] == spillover_gene$end) %>% all()
    } else {
      if_end_equals_start <- (ws286_find_spillovers$start[i] == spillover_gene$start) %>% all()
    }
    ws286_find_spillovers$should_exclude[i] <- if_end_equals_start
  } else {
    ws286_find_spillovers$should_exclude[i] <- FALSE
  }
}

ws286_find_spillovers$should_exclude %>% which() %>% length() # same number?

ws286_find_spillovers$should_exclude[is.na(ws286_find_spillovers$should_exclude)] <- "FALSE"
ws286_find_spillovers$should_exclude %<>% as.logical()
ws286_find_spillovers$should_exclude %>% which() %>% length() # still the same number?

write_xlsx(ws286_find_spillovers, "ws286_find_spillovers_ShouldYouExclSome_250bp.xlsx", col_names = T)

