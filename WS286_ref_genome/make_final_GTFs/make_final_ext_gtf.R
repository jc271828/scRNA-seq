library(tidyverse)
library(magrittr)
library(readxl)

ws286_ext_agreed_by_2reps <- read_xlsx("./ws286_BestExtAgreedBy2reps_20230414.xlsx")

ws286_ext_agreed_by_2reps$seqname %<>% as.factor() %>% fct_relevel("I", "II", "III", "IV", "V", "X", "MtDNA")

ws286_ext_agreed_by_2reps %<>% arrange(seqname, gene_id, start)

ws286_ext_agreed_by_2reps %<>% select(colnames(.)[1:9])

colnames(ws286_ext_agreed_by_2reps) <- c("#!genebuild-version WS286")

write_tsv(ws286_ext_agreed_by_2reps, "ws286_BestExtAgreedBy2reps_20230414.gtf", col_names = T, escape = "none")

