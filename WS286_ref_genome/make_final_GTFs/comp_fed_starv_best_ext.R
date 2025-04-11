library(tidyverse)
library(magrittr)
library(readxl)
library(writexl)

# load "best extension" data
fed_half_reps <- read_xlsx("../scKB_count_matrices_code_and_outputs/Fed_BestExtHalfReps_Min0.1GainRatio.xlsx", col_names = T)
starv_half_reps <- read_xlsx("../scKB_count_matrices_code_and_outputs/Starv_BestExtHalfReps_Min0.1GainRatio.xlsx", col_names = T)

# load extended annotations
ws286 <- read_xlsx("../edit_gtf_code_and_output/ws286_new.xlsx", col_names = T)
ws286_50bp <- read_xlsx("../ext_3pUtr_code_and outputs/ws286_ext_Gene3pUtrTransExon_50bp.xlsx", col_names = T)
ws286_100bp <- read_xlsx("../ext_3pUtr_code_and outputs/ws286_ext_Gene3pUtrTransExon_100bp.xlsx", col_names = T)
ws286_150bp <- read_xlsx("../ext_3pUtr_code_and outputs/ws286_ext_Gene3pUtrTransExon_150bp.xlsx", col_names = T)
ws286_200bp <- read_xlsx("../ext_3pUtr_code_and outputs/ws286_ext_Gene3pUtrTransExon_200bp.xlsx", col_names = T)
ws286_250bp <- read_xlsx("../ext_3pUtr_code_and outputs/ws286_ext_Gene3pUtrTransExon_250bp.xlsx", col_names = T)
ws286_300bp <- read_xlsx("../ext_3pUtr_code_and outputs/ws286_ext_Gene3pUtrTransExon_300bp.xlsx", col_names = T)
ws286_400bp <- read_xlsx("../ext_3pUtr_code_and outputs/ws286_ext_Gene3pUtrTransExon_400bp.xlsx", col_names = T)
ws286_500bp <- read_xlsx("../ext_3pUtr_code_and outputs/ws286_ext_Gene3pUtrTransExon_500bp.xlsx", col_names = T)

# load stop coordinates spillover data
spillover <- read_xlsx("../edit_gtf_code_and_output/ws286_find_spillovers.xlsx", col_names = T)

# now make final version of annotation
# when I changed annotations, I changed ONLY 1 gene, at most 1 three_prime_utr, at least 1 exon (stops at where gene stops) and at least 1 transcript (stops at where gene stops) for each gene_id
# what has remained unchanged throughout all annotation version: five_prime_utr, CDS, start codon, stop codon

spillover_genes <- spillover$gene_id[spillover$StopCoordSpillOver]

fed_half_reps %<>% filter(!(gene_id %in% spillover_genes))
starv_half_reps %<>% filter(!(gene_id %in% spillover_genes))

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

# make final annotation based on "agreed by at least 5 reps out of 8" for fed
ws286_to_change_fed <- filter(ws286, gene_id %in% fed_half_reps$gene_id)
ws286_to_change_fed$start_new <- NA
ws286_to_change_fed$end_new <- NA
ws286_to_change_fed$best_ext <- NA
ws286_to_change_fed$num_reps_agreed <- NA

# this for loop can take 2 hours locally
for (i in 1:nrow(ws286_to_change_fed)) {
  ws286_to_change_fed[i, "best_ext"] <- fed_half_reps %>% filter(gene_id == ws286_to_change_fed$gene_id[i]) %>% .$best_ext %>% as.numeric()
  ws286_to_change_fed[i, "num_reps_agreed"] <- fed_half_reps %>% filter(gene_id == ws286_to_change_fed$gene_id[i]) %>% .$num_reps_agreed %>% as.numeric()
}

# row order remains the same, which is crucial for doing cbind
all.equal(ws286_to_change_fed$attribute, filter(ws286_500bp, gene_id %in% fed_half_reps$gene_id) %>% .$attribute)

# cbind
ws286_to_change_fed %<>% cbind(ws286_50bp %>% filter(gene_id %in% fed_half_reps$gene_id) %>% select(c("start", "end")),
                               ws286_100bp %>% filter(gene_id %in% fed_half_reps$gene_id) %>% select(c("start", "end")),
                               ws286_150bp %>% filter(gene_id %in% fed_half_reps$gene_id) %>% select(c("start", "end")),
                               ws286_200bp %>% filter(gene_id %in% fed_half_reps$gene_id) %>% select(c("start", "end")),
                               ws286_250bp %>% filter(gene_id %in% fed_half_reps$gene_id) %>% select(c("start", "end")),
                               ws286_300bp %>% filter(gene_id %in% fed_half_reps$gene_id) %>% select(c("start", "end")),
                               ws286_400bp %>% filter(gene_id %in% fed_half_reps$gene_id) %>% select(c("start", "end")),
                               ws286_500bp %>% filter(gene_id %in% fed_half_reps$gene_id) %>% select(c("start", "end")))

colnames(ws286_to_change_fed)[25:26] %<>% paste("50bp", sep = "_")
colnames(ws286_to_change_fed)[27:28] %<>% paste("100bp", sep = "_")
colnames(ws286_to_change_fed)[29:30] %<>% paste("150bp", sep = "_")
colnames(ws286_to_change_fed)[31:32] %<>% paste("200bp", sep = "_")
colnames(ws286_to_change_fed)[33:34] %<>% paste("250bp", sep = "_")
colnames(ws286_to_change_fed)[35:36] %<>% paste("300bp", sep = "_")
colnames(ws286_to_change_fed)[37:38] %<>% paste("400bp", sep = "_")
colnames(ws286_to_change_fed)[39:40] %<>% paste("500bp", sep = "_")

ws286_to_change_fed$best_ext %<>% paste("bp", sep = "")

# this for loop doesn't take too long to run
for (i in 1:nrow(ws286_to_change_fed)) {
  best_ext <- ws286_to_change_fed$best_ext[i]
  ws286_to_change_fed[i, "start_new"] <- ws286_to_change_fed[i, paste("start", best_ext, sep = "_")]
  ws286_to_change_fed[i, "end_new"] <- ws286_to_change_fed[i, paste("end", best_ext, sep = "_")]
}

ws286$start_new <- NA
ws286$end_new <- NA
ws286$best_ext <- NA
ws286$num_reps_agreed <- NA

ws286_best_ext_fed <- ws286 %>% filter(!(gene_id %in% fed_half_reps$gene_id)) %>% rbind(ws286_to_change_fed[, 1:24])
ws286_best_ext_fed$feature %<>% as.factor() %>% fct_relevel(c("gene", "transcript", "exon", "three_prime_utr", "CDS", "start_codon", "stop_codon", "five_prime_utr"))
ws286_best_ext_fed %<>% arrange(seqname, gene_id, feature, transcript_id, exon_id, protein_id, start)

# all annotation dfs have the same row order
all.equal(ws286_best_ext_fed$attribute, ws286_sorted$attribute)
all.equal(ws286_best_ext_fed$attribute, ws286_500bp$attribute)

# clean the final annotation (best extension agreed by at least 5 out of 8 reps)
ws286_best_ext_fed <- ws286_best_ext_fed[, c(1:5, 21:24, 6:20)]
ws286_best_ext_fed$start_new[is.na(ws286_best_ext_fed$start_new)] <- ws286_best_ext_fed$start[is.na(ws286_best_ext_fed$start_new)]
ws286_best_ext_fed$end_new[is.na(ws286_best_ext_fed$end_new)] <- ws286_best_ext_fed$end[is.na(ws286_best_ext_fed$end_new)]

write_xlsx(ws286_best_ext_fed, "ws286_best_ext_fed_more_info_corrected.xlsx", col_names = T)

# how many genes were changed?
gene_start_old <- ws286_best_ext_fed %>% filter(feature == "gene") %>% .$start
gene_start_new <- ws286_best_ext_fed %>% filter(feature == "gene") %>% .$start_new
gene_end_old <- ws286_best_ext_fed %>% filter(feature == "gene") %>% .$end
gene_end_new <- ws286_best_ext_fed %>% filter(feature == "gene") %>% .$end_new

# 382 genes were changed (22 genes were extended until spilling into adjacent genes)
which(!(gene_start_old == gene_start_new)) %>% length() # 214
which(!(gene_end_old == gene_end_new)) %>% length() # 168

# 16 genes were extended by 50 bp
which(gene_start_old - 50 == gene_start_new) %>% length() # 6
which(gene_end_old + 50 == gene_end_new) %>% length() # 10

# 48 genes were extended by 100 bp
which(gene_start_old - 100 == gene_start_new) %>% length() # 28
which(gene_end_old + 100 == gene_end_new) %>% length() # 20

# 44 genes were extended by 150 bp
which(gene_start_old - 150 == gene_start_new) %>% length() # 24
which(gene_end_old + 150 == gene_end_new) %>% length() # 20

# 41 genes were extended by 200 bp
which(gene_start_old - 200 == gene_start_new) %>% length() # 18
which(gene_end_old + 200 == gene_end_new) %>% length() # 23

# 26 genes were extended by 250 bp
which(gene_start_old - 250 == gene_start_new) %>% length() # 14
which(gene_end_old + 250 == gene_end_new) %>% length() # 12

# 21 genes were extended by 300 bp
which(gene_start_old - 300 == gene_start_new) %>% length() # 14
which(gene_end_old + 300 == gene_end_new) %>% length() # 7

# 68 genes were extended by 400 bp
which(gene_start_old - 400 == gene_start_new) %>% length() # 47
which(gene_end_old + 400 == gene_end_new) %>% length() # 21

# 96 genes were extended by 500 bp
which(gene_start_old - 500 == gene_start_new) %>% length() # 53
which(gene_end_old + 500 == gene_end_new) %>% length() # 43

ws286_best_ext_fed %<>% select(-c("start", "end", "best_ext", "num_reps_agreed"))
colnames(ws286_best_ext_fed)[4:5] <- c("start", "end")

write_xlsx(ws286_best_ext_fed, "ws286_best_ext_fed_corrected.xlsx", col_names = T)

################################################################################################################################################################################################

# make final annotation based on "agreed by at least 5 reps out of 8" for starv
# now use the file ws286_sorted that was saved a bit ago
ws286 <- ws286_sorted

ws286_to_change_starv <- filter(ws286, gene_id %in% starv_half_reps$gene_id)
ws286_to_change_starv$start_new <- NA
ws286_to_change_starv$end_new <- NA
ws286_to_change_starv$best_ext <- NA
ws286_to_change_starv$num_reps_agreed <- NA

# this for loop can take 2 hours locally
for (i in 1:nrow(ws286_to_change_starv)) {
  ws286_to_change_starv[i, "best_ext"] <- starv_half_reps %>% filter(gene_id == ws286_to_change_starv$gene_id[i]) %>% .$best_ext %>% as.numeric()
  ws286_to_change_starv[i, "num_reps_agreed"] <- starv_half_reps %>% filter(gene_id == ws286_to_change_starv$gene_id[i]) %>% .$num_reps_agreed %>% as.numeric()
}

# row order remains the same, which is crucial for doing cbind
all.equal(ws286_to_change_starv$attribute, filter(ws286_sorted, gene_id %in% starv_half_reps$gene_id) %>% .$attribute)
all.equal(ws286_to_change_starv$attribute, filter(ws286_500bp, gene_id %in% starv_half_reps$gene_id) %>% .$attribute)

# cbind
ws286_to_change_starv %<>% cbind(ws286_50bp %>% filter(gene_id %in% starv_half_reps$gene_id) %>% select(c("start", "end")),
                                 ws286_100bp %>% filter(gene_id %in% starv_half_reps$gene_id) %>% select(c("start", "end")),
                                 ws286_150bp %>% filter(gene_id %in% starv_half_reps$gene_id) %>% select(c("start", "end")),
                                 ws286_200bp %>% filter(gene_id %in% starv_half_reps$gene_id) %>% select(c("start", "end")),
                                 ws286_250bp %>% filter(gene_id %in% starv_half_reps$gene_id) %>% select(c("start", "end")),
                                 ws286_300bp %>% filter(gene_id %in% starv_half_reps$gene_id) %>% select(c("start", "end")),
                                 ws286_400bp %>% filter(gene_id %in% starv_half_reps$gene_id) %>% select(c("start", "end")),
                                 ws286_500bp %>% filter(gene_id %in% starv_half_reps$gene_id) %>% select(c("start", "end")))

colnames(ws286_to_change_starv)[25:26] %<>% paste("50bp", sep = "_")
colnames(ws286_to_change_starv)[27:28] %<>% paste("100bp", sep = "_")
colnames(ws286_to_change_starv)[29:30] %<>% paste("150bp", sep = "_")
colnames(ws286_to_change_starv)[31:32] %<>% paste("200bp", sep = "_")
colnames(ws286_to_change_starv)[33:34] %<>% paste("250bp", sep = "_")
colnames(ws286_to_change_starv)[35:36] %<>% paste("300bp", sep = "_")
colnames(ws286_to_change_starv)[37:38] %<>% paste("400bp", sep = "_")
colnames(ws286_to_change_starv)[39:40] %<>% paste("500bp", sep = "_")

ws286_to_change_starv$best_ext %<>% paste("bp", sep = "")

# this for loop doesn't take too long to run
for (i in 1:nrow(ws286_to_change_starv)) {
  best_ext <- ws286_to_change_starv$best_ext[i]
  ws286_to_change_starv[i, "start_new"] <- ws286_to_change_starv[i, paste("start", best_ext, sep = "_")]
  ws286_to_change_starv[i, "end_new"] <- ws286_to_change_starv[i, paste("end", best_ext, sep = "_")]
}

ws286$start_new <- NA
ws286$end_new <- NA
ws286$best_ext <- NA
ws286$num_reps_agreed <- NA

ws286_best_ext_starv <- ws286 %>% filter(!(gene_id %in% starv_half_reps$gene_id)) %>% rbind(ws286_to_change_starv[, 1:24])
ws286_best_ext_starv$feature %<>% as.factor() %>% fct_relevel(c("gene", "transcript", "exon", "three_prime_utr", "CDS", "start_codon", "stop_codon", "five_prime_utr"))
ws286_best_ext_starv %<>% arrange(seqname, gene_id, feature, transcript_id, exon_id, protein_id, start)

# all annotation dfs have the same row order
all.equal(ws286_best_ext_starv$attribute, ws286_sorted$attribute)
all.equal(ws286_best_ext_starv$attribute, ws286_500bp$attribute)

# clean the final annotation (best extension agreed by at least 5 out of 8 reps)
ws286_best_ext_starv <- ws286_best_ext_starv[, c(1:5, 21:24, 6:20)]
ws286_best_ext_starv$start_new[is.na(ws286_best_ext_starv$start_new)] <- ws286_best_ext_starv$start[is.na(ws286_best_ext_starv$start_new)]
ws286_best_ext_starv$end_new[is.na(ws286_best_ext_starv$end_new)] <- ws286_best_ext_starv$end[is.na(ws286_best_ext_starv$end_new)]

write_xlsx(ws286_best_ext_starv, "ws286_best_ext_starv_more_info_corrected.xlsx", col_names = T)

# how many genes were changed?
gene_start_old <- ws286_best_ext_starv %>% filter(feature == "gene") %>% .$start
gene_start_new <- ws286_best_ext_starv %>% filter(feature == "gene") %>% .$start_new
gene_end_old <- ws286_best_ext_starv %>% filter(feature == "gene") %>% .$end
gene_end_new <- ws286_best_ext_starv %>% filter(feature == "gene") %>% .$end_new

# 483 genes were changed (25 genes were extended until spilling into adjacent genes)
which(!(gene_start_old == gene_start_new)) %>% length() # 253
which(!(gene_end_old == gene_end_new)) %>% length() # 230

# 16 genes were extended by 50 bp
which(gene_start_old - 50 == gene_start_new) %>% length() # 6
which(gene_end_old + 50 == gene_end_new) %>% length() # 10

# 63 genes were extended by 100 bp
which(gene_start_old - 100 == gene_start_new) %>% length() # 35
which(gene_end_old + 100 == gene_end_new) %>% length() # 28

# 58 genes were extended by 150 bp
which(gene_start_old - 150 == gene_start_new) %>% length() # 30
which(gene_end_old + 150 == gene_end_new) %>% length() # 28

# 48 genes were extended by 200 bp
which(gene_start_old - 200 == gene_start_new) %>% length() # 22
which(gene_end_old + 200 == gene_end_new) %>% length() # 26

# 25 genes were extended by 250 bp
which(gene_start_old - 250 == gene_start_new) %>% length() # 15
which(gene_end_old + 250 == gene_end_new) %>% length() # 10

# 30 genes were extended by 300 bp
which(gene_start_old - 300 == gene_start_new) %>% length() # 16
which(gene_end_old + 300 == gene_end_new) %>% length() # 14

# 84 genes were extended by 400 bp
which(gene_start_old - 400 == gene_start_new) %>% length() # 51
which(gene_end_old + 400 == gene_end_new) %>% length() # 33

# 134 genes were extended by 500 bp
which(gene_start_old - 500 == gene_start_new) %>% length() # 72
which(gene_end_old + 500 == gene_end_new) %>% length() # 62

ws286_best_ext_starv %<>% select(-c("start", "end", "best_ext", "num_reps_agreed"))
colnames(ws286_best_ext_starv)[4:5] <- c("start", "end")

write_xlsx(ws286_best_ext_starv, "ws286_best_ext_starv_corrected.xlsx", col_names = T)

# compare fed starv best extensions

# ws286_best_ext_fed <- read_xlsx("ws286_best_ext_fed_corrected.xlsx", col_names = T)
# ws286_best_ext_starv <- read_xlsx("ws286_best_ext_starv_corrected.xlsx", col_names = T)

# again this shows it's so important to have saved the sorted ws286 as ws286_sorted
ws286 <- ws286_sorted

# so it's good to do cbind
all.equal(ws286$attribute, ws286_best_ext_starv$attribute)
all.equal(ws286$attribute, ws286_best_ext_fed$attribute)

colnames(ws286)[4:5] <- c("start_raw", "end_raw")
colnames(ws286_best_ext_fed)[4:5] <- c("start_ext_fed", "end_ext_fed")
colnames(ws286_best_ext_starv)[4:5] <- c("start_ext_starv", "end_ext_starv")

comp_fed_starv <- cbind(ws286, ws286_best_ext_fed[, 4:5], ws286_best_ext_starv[, 4:5])
comp_fed_starv %<>% filter(feature == "gene") %>% select(-c("feature"))
comp_fed_starv$fed_best_ext <- NA
comp_fed_starv$starv_best_ext <- NA

for (i in 1:nrow(comp_fed_starv)) {
  comp_fed_starv[i, "fed_best_ext"] <- max(comp_fed_starv$start_raw[i] - comp_fed_starv$start_ext_fed[i], comp_fed_starv$end_ext_fed[i] - comp_fed_starv$end_raw[i])
  comp_fed_starv[i, "starv_best_ext"] <- max(comp_fed_starv$start_raw[i] - comp_fed_starv$start_ext_starv[i], comp_fed_starv$end_ext_starv[i] - comp_fed_starv$end_raw[i])
}

ggplot()+
  stat_function(color = "red", linewidth = 1, fun = function(x) x)+
  geom_point(data = comp_fed_starv, mapping = aes(x = fed_best_ext, y = starv_best_ext)) +
  theme_classic()+
  theme(aspect.ratio = 1) +
  labs(title = "best extension determined by >5/8 reps comparing fed and starv")

only_ext_in_starv <- comp_fed_starv %>% filter(fed_best_ext == 0 & starv_best_ext != 0)
only_ext_in_fed <- comp_fed_starv %>% filter(fed_best_ext != 0 & starv_best_ext == 0)

write_xlsx(only_ext_in_fed, "only_ext_in_fed.xlsx", col_names = T)
write_xlsx(only_ext_in_starv, "only_ext_in_starv.xlsx", col_names = T)




