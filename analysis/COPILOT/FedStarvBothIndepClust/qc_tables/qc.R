library(tidyverse)
library(magrittr)
library(readxl)

qc_per_sample <- read_tsv("qc_per_sample_20231006.txt")

qc_per_sample <- qc_per_sample[, c(1,2,4,8)]
colnames(qc_per_sample)[2:4] <- c("number of cells", "median number of transcripts from protein-coding genes", "number of protein-coding genes detected")
qc_per_sample$`number of cells` %<>% as.integer()
qc_per_sample$`median number of transcripts from protein-coding genes` %<>% as.integer()
qc_per_sample$`number of protein-coding genes detected` %<>% as.integer()

writexl::write_xlsx(qc_per_sample, "qc_per_sample.xlsx")

qc_per_condition <- read_tsv("qc_per_condition_20231006.txt")
