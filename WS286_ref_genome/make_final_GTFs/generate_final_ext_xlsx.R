library(tidyverse)
library(magrittr)
library(readxl)
library(writexl)

# load "best extension" data
agreed_by_2reps <- read_xlsx("../scKB_count_matrices_code_and_outputs/BestExt2Reps0.05GainRatio_20230422.xlsx", col_names = T)

# load extended annotations
ws286 <- read_xlsx("../edit_gtf_code_and_output/ws286_new.xlsx", col_names = T)
ws286_50bp <- read_xlsx("../ext_3pUtr_code_and outputs/ws286_ext_Gene3pUtrTrans_50bp_20230422.xlsx", col_names = T)
ws286_100bp <- read_xlsx("../ext_3pUtr_code_and outputs/ws286_ext_Gene3pUtrTrans_100bp_20230422.xlsx", col_names = T)
ws286_150bp <- read_xlsx("../ext_3pUtr_code_and outputs/ws286_ext_Gene3pUtrTrans_150bp_20230422.xlsx", col_names = T)
ws286_200bp <- read_xlsx("../ext_3pUtr_code_and outputs/ws286_ext_Gene3pUtrTrans_200bp_20230422.xlsx", col_names = T)
ws286_250bp <- read_xlsx("../ext_3pUtr_code_and outputs/ws286_ext_Gene3pUtrTrans_250bp_20230422.xlsx", col_names = T)
ws286_300bp <- read_xlsx("../ext_3pUtr_code_and outputs/ws286_ext_Gene3pUtrTrans_300bp_20230422.xlsx", col_names = T)
ws286_400bp <- read_xlsx("../ext_3pUtr_code_and outputs/ws286_ext_Gene3pUtrTrans_400bp_20230422.xlsx", col_names = T)
ws286_500bp <- read_xlsx("../ext_3pUtr_code_and outputs/ws286_ext_Gene3pUtrTrans_500bp_20230422.xlsx", col_names = T)

# now make final version of annotation
# when I changed annotations, I changed ONLY 1 gene, at most 1 three_prime_utr, at least 1 exon (stops at where gene stops) and at least 1 transcript (stops at where gene stops) for each gene_id
# what has remained unchanged throughout all annotation version: five_prime_utr, CDS, start codon, stop codon

arr <- function(df) {
  df$feature %<>% as.factor() %>% fct_relevel(c("gene", "transcript", "exon", "three_prime_utr", "CDS", "start_codon", "stop_codon", "five_prime_utr"))
  df %<>% arrange(seqname, gene_id, feature, transcript_id, exon_id, protein_id, start)
  return(df)
}

ws286 %<>% arr()
ws286_50bp %<>% arr()
ws286_100bp %<>% arr()
ws286_150bp %<>% arr()
ws286_200bp %<>% arr()
ws286_250bp %<>% arr()
ws286_300bp %<>% arr()
ws286_400bp %<>% arr()
ws286_500bp %<>% arr()

# all annotation dfs have the same row order 
all.equal(ws286$attribute, ws286_500bp$attribute)

# keep a copy of the sorted ws286
# because you will see later how important it is to keep the order!!!
ws286_sorted <- ws286

# make final annotation based on "agreed by at least 2 reps out of 8"
ws286_to_change_2reps <- filter(ws286, gene_id %in% agreed_by_2reps$gene_id)
ws286_to_change_2reps$start_new <- NA
ws286_to_change_2reps$end_new <- NA
ws286_to_change_2reps$best_ext <- NA
ws286_to_change_2reps$num_reps_agreed <- NA

# this for loop can take 2 hours locally
for (i in 1:nrow(ws286_to_change_2reps)) {
  ws286_to_change_2reps[i, "best_ext"] <- agreed_by_2reps %>% filter(gene_id == ws286_to_change_2reps$gene_id[i]) %>% .$best_ext %>% as.numeric()
  ws286_to_change_2reps[i, "num_reps_agreed"] <- agreed_by_2reps %>% filter(gene_id == ws286_to_change_2reps$gene_id[i]) %>% .$num_reps_agreed %>% as.numeric()
}

# row order remains the same, which is crucial for doing cbind
all.equal(ws286_to_change_2reps$attribute, filter(ws286_500bp, gene_id %in% agreed_by_2reps$gene_id) %>% .$attribute)

# cbind
ws286_to_change_2reps %<>% cbind(ws286_50bp %>% filter(gene_id %in% agreed_by_2reps$gene_id) %>% select(c("start", "end")),
                                 ws286_100bp %>% filter(gene_id %in% agreed_by_2reps$gene_id) %>% select(c("start", "end")),
                                 ws286_150bp %>% filter(gene_id %in% agreed_by_2reps$gene_id) %>% select(c("start", "end")),
                                 ws286_200bp %>% filter(gene_id %in% agreed_by_2reps$gene_id) %>% select(c("start", "end")),
                                 ws286_250bp %>% filter(gene_id %in% agreed_by_2reps$gene_id) %>% select(c("start", "end")),
                                 ws286_300bp %>% filter(gene_id %in% agreed_by_2reps$gene_id) %>% select(c("start", "end")),
                                 ws286_400bp %>% filter(gene_id %in% agreed_by_2reps$gene_id) %>% select(c("start", "end")),
                                 ws286_500bp %>% filter(gene_id %in% agreed_by_2reps$gene_id) %>% select(c("start", "end")))

colnames(ws286_to_change_2reps)[25:26] %<>% paste("50bp", sep = "_")
colnames(ws286_to_change_2reps)[27:28] %<>% paste("100bp", sep = "_")
colnames(ws286_to_change_2reps)[29:30] %<>% paste("150bp", sep = "_")
colnames(ws286_to_change_2reps)[31:32] %<>% paste("200bp", sep = "_")
colnames(ws286_to_change_2reps)[33:34] %<>% paste("250bp", sep = "_")
colnames(ws286_to_change_2reps)[35:36] %<>% paste("300bp", sep = "_")
colnames(ws286_to_change_2reps)[37:38] %<>% paste("400bp", sep = "_")
colnames(ws286_to_change_2reps)[39:40] %<>% paste("500bp", sep = "_")

ws286_to_change_2reps$best_ext %<>% paste("bp", sep = "")

# this for loop doesn't take too long to run
for (i in 1:nrow(ws286_to_change_2reps)) {
  best_ext <- ws286_to_change_2reps$best_ext[i]
  ws286_to_change_2reps[i, "start_new"] <- ws286_to_change_2reps[i, paste("start", best_ext, sep = "_")]
  ws286_to_change_2reps[i, "end_new"] <- ws286_to_change_2reps[i, paste("end", best_ext, sep = "_")]
}

# prep ws286 to generate final xlsx
ws286$start_new <- NA
ws286$end_new <- NA
ws286$best_ext <- NA
ws286$num_reps_agreed <- NA

ws286_best_ext_agreed_by_2reps <- ws286 %>% filter(!(gene_id %in% agreed_by_2reps$gene_id)) %>% rbind(ws286_to_change_2reps[, 1:24])
ws286_best_ext_agreed_by_2reps$feature %<>% as.factor() %>% fct_relevel(c("gene", "transcript", "exon", "three_prime_utr", "CDS", "start_codon", "stop_codon", "five_prime_utr"))
ws286_best_ext_agreed_by_2reps %<>% arrange(seqname, gene_id, feature, transcript_id, exon_id, protein_id, start)

# all annotation dfs have the same row order
all.equal(ws286_best_ext_agreed_by_2reps$attribute, ws286_sorted$attribute)
all.equal(ws286_best_ext_agreed_by_2reps$attribute, ws286_500bp$attribute)

# clean the final annotation (best extension agreed by at least 2 out of 8 reps)
ws286_best_ext_agreed_by_2reps <- ws286_best_ext_agreed_by_2reps[, c(1:5, 21:24, 6:20)]
ws286_best_ext_agreed_by_2reps$start_new[is.na(ws286_best_ext_agreed_by_2reps$start_new)] <- ws286_best_ext_agreed_by_2reps$start[is.na(ws286_best_ext_agreed_by_2reps$start_new)]
ws286_best_ext_agreed_by_2reps$end_new[is.na(ws286_best_ext_agreed_by_2reps$end_new)] <- ws286_best_ext_agreed_by_2reps$end[is.na(ws286_best_ext_agreed_by_2reps$end_new)]

# how many genes were changed?
gene_start_old <- ws286_best_ext_agreed_by_2reps %>% filter(feature == "gene") %>% .$start
gene_start_new <- ws286_best_ext_agreed_by_2reps %>% filter(feature == "gene") %>% .$start_new
gene_end_old <- ws286_best_ext_agreed_by_2reps %>% filter(feature == "gene") %>% .$end
gene_end_new <- ws286_best_ext_agreed_by_2reps %>% filter(feature == "gene") %>% .$end_new

# 1536 genes were changed (116 genes were extended before spilling into adjacent genes)
# 1536 is consistent with the number I got in "../scKB_count_matrices_code_and_outputs/comb_reps_get_best_ext.R"
temp1 <- which(gene_start_old != gene_start_new) %>% length()
temp2 <- which(gene_end_old != gene_end_new) %>% length()
sum(temp1, temp2) %>% paste("genes were changed", sep = " ")

# 64 genes were extended by 50 bp. 5 genes were extended by (0, 50) bp
temp1 <- which(gene_start_old - 50 == gene_start_new) %>% length()
temp2 <- which(gene_end_old + 50 == gene_end_new) %>% length()
sum(temp1, temp2) %>% paste("genes were extended by 50 bp", sep = " ")
temp1 <- which((gene_start_old - 50 < gene_start_new) & (gene_start_old - 0 > gene_start_new)) %>% length()
temp2 <- which((gene_end_old + 50 > gene_end_new) & (gene_end_old + 0 < gene_end_new)) %>% length()
sum(temp1, temp2) %>% paste("genes were extended by (0, 50) bp", sep = " ")

# 174 genes were extended by 100 bp. 21 genes were extended by (50, 100) bp
temp1 <- which(gene_start_old - 100 == gene_start_new) %>% length()
temp2 <- which(gene_end_old + 100 == gene_end_new) %>% length()
sum(temp1, temp2) %>% paste("genes were extended by 100 bp", sep = " ")
temp1 <- which((gene_start_old - 100 < gene_start_new) & (gene_start_old - 50 > gene_start_new)) %>% length()
temp2 <- which((gene_end_old + 100 > gene_end_new) & (gene_end_old + 50 < gene_end_new)) %>% length()
sum(temp1, temp2) %>% paste("genes were extended by (50, 100) bp", sep = " ")

# 150 genes were extended by 150 bp. 16 genes were extended by (100, 150) bp
temp1 <- which(gene_start_old - 150 == gene_start_new) %>% length()
temp2 <- which(gene_end_old + 150 == gene_end_new) %>% length()
sum(temp1, temp2) %>% paste("genes were extended by 150 bp", sep = " ")
temp1 <- which((gene_start_old - 150 < gene_start_new) & (gene_start_old - 100 > gene_start_new)) %>% length()
temp2 <- which((gene_end_old + 150 > gene_end_new) & (gene_end_old + 100 < gene_end_new)) %>% length()
sum(temp1, temp2) %>% paste("genes were extended by (100, 150) bp", sep = " ")

# 131 genes were extended by 200 bp. 8 genes were extended by (150, 200) bp
temp1 <- which(gene_start_old - 200 == gene_start_new) %>% length()
temp2 <- which(gene_end_old + 200 == gene_end_new) %>% length()
sum(temp1, temp2) %>% paste("genes were extended by 200 bp", sep = " ")
temp1 <- which((gene_start_old - 200 < gene_start_new) & (gene_start_old - 150 > gene_start_new)) %>% length()
temp2 <- which((gene_end_old + 200 > gene_end_new) & (gene_end_old + 150 < gene_end_new)) %>% length()
sum(temp1, temp2) %>% paste("genes were extended by (150, 200) bp", sep = " ")

# 90 genes were extended by 250 bp. 14 genes were extended by (200, 250) bp
temp1 <- which(gene_start_old - 250 == gene_start_new) %>% length()
temp2 <- which(gene_end_old + 250 == gene_end_new) %>% length()
sum(temp1, temp2) %>% paste("genes were extended by 250 bp", sep = " ")
temp1 <- which((gene_start_old - 250 < gene_start_new) & (gene_start_old - 200 > gene_start_new)) %>% length()
temp2 <- which((gene_end_old + 250 > gene_end_new) & (gene_end_old + 200 < gene_end_new)) %>% length()
sum(temp1, temp2) %>% paste("genes were extended by (200, 250) bp", sep = " ")

# 98 genes were extended by 300 bp. 9 genes were extended by (250, 300) bp
temp1 <- which(gene_start_old - 300 == gene_start_new) %>% length()
temp2 <- which(gene_end_old + 300 == gene_end_new) %>% length()
sum(temp1, temp2) %>% paste("genes were extended by 300 bp", sep = " ")
temp1 <- which((gene_start_old - 300 < gene_start_new) & (gene_start_old - 250 > gene_start_new)) %>% length()
temp2 <- which((gene_end_old + 300 > gene_end_new) & (gene_end_old + 250 < gene_end_new)) %>% length()
sum(temp1, temp2) %>% paste("genes were extended by (250, 300) bp", sep = " ")

# 253 genes were extended by 400 bp. 25 genes were extended by (300, 400) bp
temp1 <- which(gene_start_old - 400 == gene_start_new) %>% length()
temp2 <- which(gene_end_old + 400 == gene_end_new) %>% length()
sum(temp1, temp2) %>% paste("genes were extended by 400 bp", sep = " ")
temp1 <- which((gene_start_old - 400 < gene_start_new) & (gene_start_old - 300 > gene_start_new)) %>% length()
temp2 <- which((gene_end_old + 400 > gene_end_new) & (gene_end_old + 300 < gene_end_new)) %>% length()
sum(temp1, temp2) %>% paste("genes were extended by (300, 400) bp", sep = " ")

# 460 genes were extended by 500 bp. 18 genes were extended by (400, 500) bp
temp1 <- which(gene_start_old - 500 == gene_start_new) %>% length()
temp2 <- which(gene_end_old + 500 == gene_end_new) %>% length()
sum(temp1, temp2) %>% paste("genes were extended by 500 bp", sep = " ")
temp1 <- which((gene_start_old - 500 < gene_start_new) & (gene_start_old - 400 > gene_start_new)) %>% length()
temp2 <- which((gene_end_old + 500 > gene_end_new) & (gene_end_old + 400 < gene_end_new)) %>% length()
sum(temp1, temp2) %>% paste("genes were extended by (400, 500) bp", sep = " ")

write_xlsx(ws286_best_ext_agreed_by_2reps, "ws286_BestExtAgreedBy2reps_more_info_20230422.xlsx", col_names = T)

ws286_best_ext_agreed_by_2reps %<>% select(-c("start", "end", "best_ext", "num_reps_agreed"))
colnames(ws286_best_ext_agreed_by_2reps)[4:5] <- c("start", "end")

write_xlsx(ws286_best_ext_agreed_by_2reps, "ws286_BestExtAgreedBy2reps_20230422.xlsx", col_names = T)


















