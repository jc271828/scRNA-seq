library(tidyverse)
library(magrittr)
library(readxl)

setwd("after_rlzing_DCC_does_weird_sorting/")

SO_raw <- read_xlsx("../../edit_gtf_code_and_output/ws286_find_spillovers_ShouldYouExclSome.xlsx")
SO_50bp <- read_xlsx("ws286_find_spillovers_ShouldYouExclSome_50bp.xlsx")
SO_100bp <- read_xlsx("ws286_find_spillovers_ShouldYouExclSome_100bp.xlsx")
SO_150bp <- read_xlsx("ws286_find_spillovers_ShouldYouExclSome_150bp.xlsx")
SO_200bp <- read_xlsx("ws286_find_spillovers_ShouldYouExclSome_200bp.xlsx")
SO_250bp <- read_xlsx("ws286_find_spillovers_ShouldYouExclSome_250bp.xlsx")
SO_300bp <- read_xlsx("ws286_find_spillovers_ShouldYouExclSome_300bp.xlsx")
SO_400bp <- read_xlsx("ws286_find_spillovers_ShouldYouExclSome_400bp.xlsx")
SO_500bp <- read_xlsx("ws286_find_spillovers_ShouldYouExclSome_500bp.xlsx")

# define a function to find whether Stop Coordinate Spillovers acutally equal gene start/end instead of spilling over INTO other genes
if_spill_over_open_interval <- function(df) {
  df$equal_endpoint <- NA
  
  for (i in which(df$StopCoordSpillOver)) {
    SO_genes <- df$SpillIntoWhichGenes[i] %>% str_split(",") %>% unlist()
    endpoints <- filter(df, gene_id %in% SO_genes) %>% .[, c("start", "end")] %>% unlist()
    
    if (df$strand[i] == "+") {
      df$equal_endpoint[i] <- (sum(df$end[i] == endpoints) == length(SO_genes))
    } else {
      df$equal_endpoint[i] <- (sum(df$start[i] == endpoints) == length(SO_genes))
    }
  }
  
  df[is.na(df$equal_endpoint), "equal_endpoint"] <- FALSE
  
  return(df)
}

SO_raw %<>% if_spill_over_open_interval()
SO_50bp %<>% if_spill_over_open_interval()
SO_100bp %<>% if_spill_over_open_interval()
SO_150bp %<>% if_spill_over_open_interval()
SO_200bp %<>% if_spill_over_open_interval()
SO_250bp %<>% if_spill_over_open_interval()
SO_300bp %<>% if_spill_over_open_interval()
SO_400bp %<>% if_spill_over_open_interval()
SO_500bp %<>% if_spill_over_open_interval()

# are there cases where a gene spills over into multiple genes but its stop coordinate all equals endpoints?
# no
SO_raw %>% filter(more_than_one_spillover == TRUE & equal_endpoint == TRUE)
SO_50bp %>% filter(more_than_one_spillover == TRUE & equal_endpoint == TRUE)
SO_100bp %>% filter(more_than_one_spillover == TRUE & equal_endpoint == TRUE)
SO_150bp %>% filter(more_than_one_spillover == TRUE & equal_endpoint == TRUE)
SO_200bp %>% filter(more_than_one_spillover == TRUE & equal_endpoint == TRUE)
SO_250bp %>% filter(more_than_one_spillover == TRUE & equal_endpoint == TRUE)

# see if equal_endpoint >= should_exclude
# yes
sum(SO_raw$equal_endpoint)
sum(SO_raw$should_exclude)

sum(SO_50bp$equal_endpoint)
sum(SO_50bp$should_exclude)

sum(SO_100bp$equal_endpoint)
sum(SO_100bp$should_exclude)

sum(SO_150bp$equal_endpoint)
sum(SO_150bp$should_exclude)

sum(SO_200bp$equal_endpoint)
sum(SO_200bp$should_exclude)

sum(SO_250bp$equal_endpoint)
sum(SO_250bp$should_exclude)

# see if Open Interval "Real" spillover cases stay the same
sum(SO_raw$StopCoordSpillOver) - sum(SO_raw$equal_endpoint) # 18451
sum(SO_50bp$StopCoordSpillOver) - sum(SO_50bp$equal_endpoint) # 18451
sum(SO_100bp$StopCoordSpillOver) - sum(SO_100bp$equal_endpoint) # 18451
sum(SO_150bp$StopCoordSpillOver) - sum(SO_150bp$equal_endpoint) # 18451
sum(SO_200bp$StopCoordSpillOver) - sum(SO_200bp$equal_endpoint) # 18451
sum(SO_250bp$StopCoordSpillOver) - sum(SO_250bp$equal_endpoint) # 18451
sum(SO_300bp$StopCoordSpillOver) - sum(SO_300bp$equal_endpoint) # 18451
sum(SO_400bp$StopCoordSpillOver) - sum(SO_400bp$equal_endpoint) # 18451
sum(SO_500bp$StopCoordSpillOver) - sum(SO_500bp$equal_endpoint) # 18451




