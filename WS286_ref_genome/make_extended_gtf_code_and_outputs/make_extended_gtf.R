library(tidyverse)
library(magrittr)
library(readxl)

ws286 <- read_xlsx("../edit_gtf/ws286_new.xlsx", col_names = T)
ws286_050bp <- read_xlsx("../extend_3pUtr/ws286_ext_Gene_PriTrans_CorrLastExon_3pUtr_50bp.xlsx")

ws286_050bp_DidntExtExon <- read_xlsx("../extend_3pUtr/ws286_ext_Gene_PriTrans_3pUtr_50bp.xlsx")
ws286_100bp <- read_xlsx("../extend_3pUtr/ws286_ext_Gene_PriTrans_CorrLastExon_3pUtr_100bp.xlsx", col_names = T)
ws286_150bp <- read_xlsx("../extend_3pUtr/ws286_ext_Gene_PriTrans_CorrLastExon_3pUtr_150bp.xlsx", col_names = T)
ws286_200bp <- read_xlsx("../extend_3pUtr/ws286_ext_Gene_PriTrans_CorrLastExon_3pUtr_200bp.xlsx", col_names = T)
ws286_250bp <- read_xlsx("../extend_3pUtr/ws286_ext_Gene_PriTrans_CorrLastExon_3pUtr_250bp.xlsx", col_names = T)
ws286_300bp <- read_xlsx("../extend_3pUtr/ws286_ext_Gene_PriTrans_CorrLastExon_3pUtr_300bp.xlsx", col_names = T)
ws286_400bp <- read_xlsx("../extend_3pUtr/ws286_ext_Gene_PriTrans_CorrLastExon_3pUtr_400bp.xlsx", col_names = T)
ws286_500bp <- read_xlsx("../extend_3pUtr/ws286_ext_Gene_PriTrans_CorrLastExon_3pUtr_500bp.xlsx", col_names = T)

ws286$seqname %<>% as.factor() %>% fct_relevel("I", "II", "III", "IV", "V", "X", "MtDNA")
ws286_050bp$seqname %<>% as.factor() %>% fct_relevel("I", "II", "III", "IV", "V", "X", "MtDNA")

ws286_050bp_DidntExtExon$seqname %<>% as.factor() %>% fct_relevel("I", "II", "III", "IV", "V", "X", "MtDNA")
ws286_100bp$seqname %<>% as.factor() %>% fct_relevel("I", "II", "III", "IV", "V", "X", "MtDNA")
ws286_150bp$seqname %<>% as.factor() %>% fct_relevel("I", "II", "III", "IV", "V", "X", "MtDNA")
ws286_200bp$seqname %<>% as.factor() %>% fct_relevel("I", "II", "III", "IV", "V", "X", "MtDNA")
ws286_250bp$seqname %<>% as.factor() %>% fct_relevel("I", "II", "III", "IV", "V", "X", "MtDNA")
ws286_300bp$seqname %<>% as.factor() %>% fct_relevel("I", "II", "III", "IV", "V", "X", "MtDNA")
ws286_400bp$seqname %<>% as.factor() %>% fct_relevel("I", "II", "III", "IV", "V", "X", "MtDNA")
ws286_500bp$seqname %<>% as.factor() %>% fct_relevel("I", "II", "III", "IV", "V", "X", "MtDNA")


ws286 %<>% arrange(seqname, gene_id, gene_name, transcript_id, start, end, feature)
ws286_050bp %<>% arrange(seqname, gene_id, gene_name, transcript_id, start, end, feature)

ws286_050bp_DidntExtExon %<>% arrange(seqname, gene_id, gene_name, transcript_id, start, end, feature)
ws286_100bp %<>% arrange(seqname, gene_id, gene_name, transcript_id, start, end, feature)
ws286_150bp %<>% arrange(seqname, gene_id, gene_name, transcript_id, start, end, feature)
ws286_200bp %<>% arrange(seqname, gene_id, gene_name, transcript_id, start, end, feature)
ws286_250bp %<>% arrange(seqname, gene_id, gene_name, transcript_id, start, end, feature)
ws286_300bp %<>% arrange(seqname, gene_id, gene_name, transcript_id, start, end, feature)
ws286_400bp %<>% arrange(seqname, gene_id, gene_name, transcript_id, start, end, feature)
ws286_500bp %<>% arrange(seqname, gene_id, gene_name, transcript_id, start, end, feature)

ws286 %<>% select(colnames(.)[1:9])
ws286_050bp %<>% select(colnames(.)[1:9])

ws286_050bp_DidntExtExon %<>% select(colnames(.)[1:9])
ws286_100bp %<>% select(colnames(.)[1:9])
ws286_150bp %<>% select(colnames(.)[1:9])
ws286_200bp %<>% select(colnames(.)[1:9])
ws286_250bp %<>% select(colnames(.)[1:9])
ws286_300bp %<>% select(colnames(.)[1:9])
ws286_400bp %<>% select(colnames(.)[1:9])
ws286_500bp %<>% select(colnames(.)[1:9])

colnames(ws286) <- c("#!genebuild-version WS286")
colnames(ws286_050bp) <- c("#!genebuild-version WS286")

colnames(ws286_050bp_DidntExtExon) <- c("#!genebuild-version WS286")
colnames(ws286_100bp) <- c("#!genebuild-version WS286")
colnames(ws286_150bp) <- c("#!genebuild-version WS286")
colnames(ws286_200bp) <- c("#!genebuild-version WS286")
colnames(ws286_250bp) <- c("#!genebuild-version WS286")
colnames(ws286_300bp) <- c("#!genebuild-version WS286")
colnames(ws286_400bp) <- c("#!genebuild-version WS286")
colnames(ws286_500bp) <- c("#!genebuild-version WS286")

write_tsv(ws286, "ws286_sorted.gtf", col_names = T, escape = "none")
write_tsv(ws286_050bp, "ws286_ext_Gene_PriTrans_CorrLastExon_3pUtr_50bp.gtf", col_names = T, escape = "none")

write_tsv(ws286_050bp_DidntExtExon, "ws286_ext_Gene_PriTrans_3pUtr_50bp.gtf", col_names = T, escape = "none")
write_tsv(ws286_100bp, "ws286_ext_Gene_PriTrans_CorrLastExon_3pUtr_100bp.gtf", col_names = T, escape = "none")
write_tsv(ws286_150bp, "ws286_ext_Gene_PriTrans_CorrLastExon_3pUtr_150bp.gtf", col_names = T, escape = "none")
write_tsv(ws286_200bp, "ws286_ext_Gene_PriTrans_CorrLastExon_3pUtr_200bp.gtf", col_names = T, escape = "none")
write_tsv(ws286_250bp, "ws286_ext_Gene_PriTrans_CorrLastExon_3pUtr_250bp.gtf", col_names = T, escape = "none")
write_tsv(ws286_300bp, "ws286_ext_Gene_PriTrans_CorrLastExon_3pUtr_300bp.gtf", col_names = T, escape = "none")
write_tsv(ws286_400bp, "ws286_ext_Gene_PriTrans_CorrLastExon_3pUtr_400bp.gtf", col_names = T, escape = "none")
write_tsv(ws286_500bp, "ws286_ext_Gene_PriTrans_CorrLastExon_3pUtr_500bp.gtf", col_names = T, escape = "none")



