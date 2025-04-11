# takes 30m on the cluster with 8 cpus and 48G mem
# consumes a LOT of memory, don't run this locally
library(tidyverse)
library(magrittr)
library(readxl)
library(writexl)

ws260 <- read_xlsx("ws260.xlsx", col_names = T)
ws260$gene_id <- NA
ws260$gene_version <- NA
ws260$gene_source <- NA
ws260$gene_biotype <- NA
ws260$gene_name <- NA
ws260$transcript_id <- NA
ws260$transcript_source <- NA
ws260$transcript_biotype <- NA
ws260$protein_id <- NA
ws260$exon_id <- NA
ws260$exon_number <- NA


replace_non_na <- function(x, string) {
  string[str_detect(string, x)] %>% ifelse(length(.) == 0, NA, .)
}

DfSplit <- function(y) {
  temp <- ws260[(y-9999): y, ]
  
  for (i in 1:10000) {
    string <- str_split(temp$attribute[i], ";") %>% unlist()
    
    temp$gene_id[i] <- replace_non_na("gene_id", string)
    temp$gene_version[i] <- replace_non_na("gene_version", string)
    temp$gene_source[i] <- replace_non_na("gene_source", string)
    temp$gene_biotype[i] <- replace_non_na("gene_biotype", string)
    temp$gene_name[i] <- replace_non_na("gene_name", string)
    temp$transcript_id[i] <- replace_non_na("transcript_id", string)
    temp$transcript_source[i] <- replace_non_na("transcript_source", string)
    temp$transcript_biotype[i] <- replace_non_na("transcript_biotype", string)
    temp$protein_id[i] <- replace_non_na("protein_id", string)
    temp$exon_id[i] <- replace_non_na("exon_id", string)
    temp$exon_number[i] <- replace_non_na("exon_number", string)
  }
  return(temp)
}

nrow(ws260)
ws260_new <- NA

for (j in 1:72) {
  ws260_new %<>% rbind(DfSplit(10000 * j))
}

ws260_new <- ws260_new[-1, ]

temp <- ws260[720001: nrow(ws260), ]
for (i in 1:(nrow(ws260) - 720000)) {
  string <- str_split(temp$attribute[i], ";") %>% unlist()
  
  temp$gene_id[i] <- replace_non_na("gene_id", string)
  temp$gene_version[i] <- replace_non_na("gene_version", string)
  temp$gene_source[i] <- replace_non_na("gene_source", string)
  temp$gene_biotype[i] <- replace_non_na("gene_biotype", string)
  temp$gene_name[i] <- replace_non_na("gene_name", string)
  temp$transcript_id[i] <- replace_non_na("transcript_id", string)
  temp$transcript_source[i] <- replace_non_na("transcript_source", string)
  temp$transcript_biotype[i] <- replace_non_na("transcript_biotype", string)
  temp$protein_id[i] <- replace_non_na("protein_id", string)
  temp$exon_id[i] <- replace_non_na("exon_id", string)
  temp$exon_number[i] <- replace_non_na("exon_number", string)
}

ws260_new %<>% rbind(temp)
rm(list = c("temp", "ws260"))

# ws260_new <- read_xlsx("ws260_new.xlsx", col_names = T)

ws260_new$gene_id %<>% str_remove_all("gene_id |\"| ")
ws260_new$gene_version %<>% str_remove_all("gene_version |\"| ")
ws260_new$gene_source %<>% str_remove_all("gene_source |\"| ")
ws260_new$gene_biotype %<>% str_remove_all("gene_biotype |\"| ")
ws260_new$gene_name %<>% str_remove_all("gene_name |\"| ")
ws260_new$transcript_id %<>% str_remove_all("transcript_id |\"| ")
ws260_new$transcript_source %<>% str_remove_all("transcript_source |\"| ")
ws260_new$transcript_biotype %<>% str_remove_all("transcript_biotype |\"| ")
ws260_new$protein_id %<>% str_remove_all("protein_id |\"| ")
ws260_new$exon_id %<>% str_remove_all("exon_id |\"| ")
ws260_new$exon_number %<>% str_remove_all("exon_number |\"| ")

# overwrites the current ws826_new.xlsx in the WD
write_xlsx(ws260_new, "ws260_new.xlsx", col_names = T)
