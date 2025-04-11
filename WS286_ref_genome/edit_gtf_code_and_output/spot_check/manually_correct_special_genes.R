library(tidyverse)
library(magrittr)
library(readxl)
library(writexl)

spillover <- read_xlsx("ws286_find_spillovers_ShouldYouExclSome.xlsx")

spillover$adjacent_pos <- NA
spillover$adjacent_strand <- NA
spillover$adjacent_gene <- NA

for (i in which(spillover$should_exclude)) {
  if (spillover$strand[i] == "+") {
    next_gene <- filter(spillover, seqname == spillover$seqname[i] & start > spillover$end[i]) %>% arrange(start, desc(strand)) %>% .[1, ]
    if (nrow(next_gene) > 0) {
      spillover$adjacent_pos[i] <- next_gene$start
      spillover$adjacent_strand[i] <- next_gene$strand
      spillover$adjacent_gene[i] <- next_gene$gene_id
    } else {
      spillover$adjacent_pos[i] <- NA
      spillover$adjacent_strand[i] <- NA
      spillover$adjacent_gene[i] <- NA
    }
  } else {
    previous_gene <- filter(spillover, seqname == spillover$seqname[i] & end < spillover$start[i]) %>% arrange(desc(end), strand) %>% .[1, ]
    if (nrow(previous_gene) > 0) {
      spillover$adjacent_pos[i] <- previous_gene$end
      spillover$adjacent_strand[i] <- previous_gene$strand
      spillover$adjacent_gene[i] <- previous_gene$gene_id
    } else {
      spillover$adjacent_pos[i] <- NA
      spillover$adjacent_strand[i] <- NA
      spillover$adjacent_gene[i] <- NA
    }
  }
}

temp <- spillover %>% filter(should_exclude == TRUE)
temp %<>% select("seqname", "start", "end", "strand", "gene_id", "StopCoordSpillOver", "SpillIntoWhichGenes", "more_than_one_spillover", "should_exclude", "adjacent_pos", "adjacent_strand", "adjacent_gene")

temp$dis <- NA
for (i in 1:nrow(temp)) {
  if (temp$strand[i] == "+") {
    temp$dis[i] <- temp$adjacent_pos[i] - temp$end[i]
  } else {
    temp$dis[i] <- temp$start[i] - temp$adjacent_pos[i]
  }
}

temp %<>% filter(dis < 1000)
temp$same_strand <- NA
temp$same_strand <- (temp$strand == temp$adjacent_strand)
temp %<>% filter(same_strand == FALSE)

temp %<>% select(-c("dis", "same_strand"))

# temp %>% filter(adjacent_gene %in% c("WBGene00199432", "WBGene00196525", "WBGene00196241"))

temp$adjacent_50bp <- NA
temp$adjacent_100bp <- NA
temp$adjacent_150bp <- NA
temp$adjacent_200bp <- NA
temp$adjacent_250bp <- NA
temp$adjacent_300bp <- NA
temp$adjacent_400bp <- NA
temp$adjacent_500bp <- NA

temp$gene_50bp <- NA
temp$gene_100bp <- NA
temp$gene_150bp <- NA
temp$gene_200bp <- NA
temp$gene_250bp <- NA
temp$gene_300bp <- NA
temp$gene_400bp <- NA
temp$gene_500bp <- NA

temp$gene_pos <- NA
for (i in 1:nrow(temp)) {
  if(temp$strand[i] == "+") {
    temp$gene_pos[i] <- temp$end[i]
  } else {
    temp$gene_pos[i] <- temp$start[i]
  }
}

temp %<>% select(-c("start", "end"))

# ext to 50bp
ext <- 50
for (i in 1:nrow(temp)) {
  if(temp$adjacent_strand[i] == "+") {
    if (temp$adjacent_pos[i] + (2*ext) < temp$gene_pos[i]) {
      temp$adjacent_50bp[i] <- temp$adjacent_pos[i] + ext
      temp$gene_50bp[i] <- temp$gene_pos[i] - ext
    } else {
      dis <- temp$gene_pos[i] - temp$adjacent_pos[i]
      temp$adjacent_50bp[i] <- temp$adjacent_pos[i] + floor(dis/2)
      temp$gene_50bp[i] <- temp$gene_pos[i] - floor(dis/2)
    }
  } else {
    if (temp$adjacent_pos[i] - (2*ext) > temp$gene_pos[i]) {
      temp$adjacent_50bp[i] <- temp$adjacent_pos[i] - ext
      temp$gene_50bp[i] <- temp$gene_pos[i] + ext
    } else {
      dis <- temp$adjacent_pos[i] - temp$gene_pos[i]
      temp$adjacent_50bp[i] <- temp$adjacent_pos[i] - floor(dis/2)
      temp$gene_50bp[i] <- temp$gene_pos[i] + floor(dis/2)
    }
  }
}

# ext to 100bp
ext <- 50
for (i in 1:nrow(temp)) {
  if(temp$adjacent_strand[i] == "+") {
    if (temp$adjacent_50bp[i] + (2*ext) < temp$gene_50bp[i]) {
      temp$adjacent_100bp[i] <- temp$adjacent_50bp[i] + ext
      temp$gene_100bp[i] <- temp$gene_50bp[i] - ext
    } else {
      dis <- temp$gene_50bp[i] - temp$adjacent_50bp[i]
      temp$adjacent_100bp[i] <- temp$adjacent_50bp[i] + floor(dis/2)
      temp$gene_100bp[i] <- temp$gene_50bp[i] - floor(dis/2)
    }
  } else {
    if (temp$adjacent_50bp[i] - (2*ext) > temp$gene_50bp[i]) {
      temp$adjacent_100bp[i] <- temp$adjacent_50bp[i] - ext
      temp$gene_100bp[i] <- temp$gene_50bp[i] + ext
    } else {
      dis <- temp$adjacent_50bp[i] - temp$gene_50bp[i]
      temp$adjacent_100bp[i] <- temp$adjacent_50bp[i] - floor(dis/2)
      temp$gene_100bp[i] <- temp$gene_50bp[i] + floor(dis/2)
    }
  }
}

# ext to 150bp
ext <- 50
for (i in 1:nrow(temp)) {
  if(temp$adjacent_strand[i] == "+") {
    if (temp$adjacent_100bp[i] + (2*ext) < temp$gene_100bp[i]) {
      temp$adjacent_150bp[i] <- temp$adjacent_100bp[i] + ext
      temp$gene_150bp[i] <- temp$gene_100bp[i] - ext
    } else {
      dis <- temp$gene_100bp[i] - temp$adjacent_100bp[i]
      temp$adjacent_150bp[i] <- temp$adjacent_100bp[i] + floor(dis/2)
      temp$gene_150bp[i] <- temp$gene_100bp[i] - floor(dis/2)
    }
  } else {
    if (temp$adjacent_100bp[i] - (2*ext) > temp$gene_100bp[i]) {
      temp$adjacent_150bp[i] <- temp$adjacent_100bp[i] - ext
      temp$gene_150bp[i] <- temp$gene_100bp[i] + ext
    } else {
      dis <- temp$adjacent_100bp[i] - temp$gene_100bp[i]
      temp$adjacent_150bp[i] <- temp$adjacent_100bp[i] - floor(dis/2)
      temp$gene_150bp[i] <- temp$gene_100bp[i] + floor(dis/2)
    }
  }
}

# ext to 200bp
ext <- 50
for (i in 1:nrow(temp)) {
  if(temp$adjacent_strand[i] == "+") {
    if (temp$adjacent_150bp[i] + (2*ext) < temp$gene_150bp[i]) {
      temp$adjacent_200bp[i] <- temp$adjacent_150bp[i] + ext
      temp$gene_200bp[i] <- temp$gene_150bp[i] - ext
    } else {
      dis <- temp$gene_150bp[i] - temp$adjacent_150bp[i]
      temp$adjacent_200bp[i] <- temp$adjacent_150bp[i] + floor(dis/2)
      temp$gene_200bp[i] <- temp$gene_150bp[i] - floor(dis/2)
    }
  } else {
    if (temp$adjacent_150bp[i] - (2*ext) > temp$gene_150bp[i]) {
      temp$adjacent_200bp[i] <- temp$adjacent_150bp[i] - ext
      temp$gene_200bp[i] <- temp$gene_150bp[i] + ext
    } else {
      dis <- temp$adjacent_150bp[i] - temp$gene_150bp[i]
      temp$adjacent_200bp[i] <- temp$adjacent_150bp[i] - floor(dis/2)
      temp$gene_200bp[i] <- temp$gene_150bp[i] + floor(dis/2)
    }
  }
}

# ext to 250bp
ext <- 50
for (i in 1:nrow(temp)) {
  if(temp$adjacent_strand[i] == "+") {
    if (temp$adjacent_200bp[i] + (2*ext) < temp$gene_200bp[i]) {
      temp$adjacent_250bp[i] <- temp$adjacent_200bp[i] + ext
      temp$gene_250bp[i] <- temp$gene_200bp[i] - ext
    } else {
      dis <- temp$gene_200bp[i] - temp$adjacent_200bp[i]
      temp$adjacent_250bp[i] <- temp$adjacent_200bp[i] + floor(dis/2)
      temp$gene_250bp[i] <- temp$gene_200bp[i] - floor(dis/2)
    }
  } else {
    if (temp$adjacent_200bp[i] - (2*ext) > temp$gene_200bp[i]) {
      temp$adjacent_250bp[i] <- temp$adjacent_200bp[i] - ext
      temp$gene_250bp[i] <- temp$gene_200bp[i] + ext
    } else {
      dis <- temp$adjacent_200bp[i] - temp$gene_200bp[i]
      temp$adjacent_250bp[i] <- temp$adjacent_200bp[i] - floor(dis/2)
      temp$gene_250bp[i] <- temp$gene_200bp[i] + floor(dis/2)
    }
  }
}

# ext to 300bp
ext <- 50
for (i in 1:nrow(temp)) {
  if(temp$adjacent_strand[i] == "+") {
    if (temp$adjacent_250bp[i] + (2*ext) < temp$gene_250bp[i]) {
      temp$adjacent_300bp[i] <- temp$adjacent_250bp[i] + ext
      temp$gene_300bp[i] <- temp$gene_250bp[i] - ext
    } else {
      dis <- temp$gene_250bp[i] - temp$adjacent_250bp[i]
      temp$adjacent_300bp[i] <- temp$adjacent_250bp[i] + floor(dis/2)
      temp$gene_300bp[i] <- temp$gene_250bp[i] - floor(dis/2)
    }
  } else {
    if (temp$adjacent_250bp[i] - (2*ext) > temp$gene_250bp[i]) {
      temp$adjacent_300bp[i] <- temp$adjacent_250bp[i] - ext
      temp$gene_300bp[i] <- temp$gene_250bp[i] + ext
    } else {
      dis <- temp$adjacent_250bp[i] - temp$gene_250bp[i]
      temp$adjacent_300bp[i] <- temp$adjacent_250bp[i] - floor(dis/2)
      temp$gene_300bp[i] <- temp$gene_250bp[i] + floor(dis/2)
    }
  }
}

# ext to 400bp
ext <- 100
for (i in 1:nrow(temp)) {
  if(temp$adjacent_strand[i] == "+") {
    if (temp$adjacent_300bp[i] + (2*ext) < temp$gene_300bp[i]) {
      temp$adjacent_400bp[i] <- temp$adjacent_300bp[i] + ext
      temp$gene_400bp[i] <- temp$gene_300bp[i] - ext
    } else {
      dis <- temp$gene_300bp[i] - temp$adjacent_300bp[i]
      temp$adjacent_400bp[i] <- temp$adjacent_300bp[i] + floor(dis/2)
      temp$gene_400bp[i] <- temp$gene_300bp[i] - floor(dis/2)
    }
  } else {
    if (temp$adjacent_300bp[i] - (2*ext) > temp$gene_300bp[i]) {
      temp$adjacent_400bp[i] <- temp$adjacent_300bp[i] - ext
      temp$gene_400bp[i] <- temp$gene_300bp[i] + ext
    } else {
      dis <- temp$adjacent_300bp[i] - temp$gene_300bp[i]
      temp$adjacent_400bp[i] <- temp$adjacent_300bp[i] - floor(dis/2)
      temp$gene_400bp[i] <- temp$gene_300bp[i] + floor(dis/2)
    }
  }
}

# ext to 500bp
ext <- 100
for (i in 1:nrow(temp)) {
  if(temp$adjacent_strand[i] == "+") {
    if (temp$adjacent_400bp[i] + (2*ext) < temp$gene_400bp[i]) {
      temp$adjacent_500bp[i] <- temp$adjacent_400bp[i] + ext
      temp$gene_500bp[i] <- temp$gene_400bp[i] - ext
    } else {
      dis <- temp$gene_400bp[i] - temp$adjacent_400bp[i]
      temp$adjacent_500bp[i] <- temp$adjacent_400bp[i] + floor(dis/2)
      temp$gene_500bp[i] <- temp$gene_400bp[i] - floor(dis/2)
    }
  } else {
    if (temp$adjacent_400bp[i] - (2*ext) > temp$gene_400bp[i]) {
      temp$adjacent_500bp[i] <- temp$adjacent_400bp[i] - ext
      temp$gene_500bp[i] <- temp$gene_400bp[i] + ext
    } else {
      dis <- temp$adjacent_400bp[i] - temp$gene_400bp[i]
      temp$adjacent_500bp[i] <- temp$adjacent_400bp[i] - floor(dis/2)
      temp$gene_500bp[i] <- temp$gene_400bp[i] + floor(dis/2)
    }
  }
}

write_xlsx(temp, "manually_corrected_special_genes.xlsx")

hh <- temp[, c(3, 8, 10, 11:27, 2, 9)]
colnames(hh)
hh1 <- hh[, c(1, 21, 20, 12:19)]
hh2 <- hh[, c(3, 22, 2, 4:11)]
colnames(hh2) <- colnames(hh1)

hh <- rbind(hh1, hh2)

# write_xlsx(hh, "manually_corrected_special_genes_less_info.xlsx")
# hh <- read_xlsx("manually_corrected_special_genes_less_info.xlsx")

hh$start <- NA
hh$start_50bp <- NA
hh$start_100bp <- NA
hh$start_150bp <- NA
hh$start_200bp <- NA
hh$start_250bp <- NA
hh$start_300bp <- NA
hh$start_400bp <- NA
hh$start_500bp <- NA

hh$end <- NA
hh$end_50bp <- NA
hh$end_100bp <- NA
hh$end_150bp <- NA
hh$end_200bp <- NA
hh$end_250bp <- NA
hh$end_300bp <- NA
hh$end_400bp <- NA
hh$end_500bp <- NA

for (i in 1:nrow(hh)) {
  if (hh$strand[i] == "+") {
    hh$end[i] <- hh$gene_pos[i]
    hh$end_50bp[i] <- hh$gene_50bp[i]
    hh$end_100bp[i] <- hh$gene_100bp[i]
    hh$end_150bp[i] <- hh$gene_150bp[i]
    hh$end_200bp[i] <- hh$gene_200bp[i]
    hh$end_250bp[i] <- hh$gene_250bp[i]
    hh$end_300bp[i] <- hh$gene_300bp[i]
    hh$end_400bp[i] <- hh$gene_400bp[i]
    hh$end_500bp[i] <- hh$gene_500bp[i]
  } else {
    hh$start[i] <- hh$gene_pos[i]
    hh$start_50bp[i] <- hh$gene_50bp[i]
    hh$start_100bp[i] <- hh$gene_100bp[i]
    hh$start_150bp[i] <- hh$gene_150bp[i]
    hh$start_200bp[i] <- hh$gene_200bp[i]
    hh$start_250bp[i] <- hh$gene_250bp[i]
    hh$start_300bp[i] <- hh$gene_300bp[i]
    hh$start_400bp[i] <- hh$gene_400bp[i]
    hh$start_500bp[i] <- hh$gene_500bp[i]
  }
}

hh <- hh[, c(1:2, 12:29)]

setwd("../ext_3pUtr_code_and outputs/")
raw <- read_xlsx("../edit_gtf_code_and_output/ws286_new.xlsx")
ext_50bp <- read_xlsx("ws286_ext_Gene3pUtrTrans_50bp_20230422.xlsx")
ext_100bp <- read_xlsx("ws286_ext_Gene3pUtrTrans_100bp_20230422.xlsx")
ext_150bp <- read_xlsx("ws286_ext_Gene3pUtrTrans_150bp_20230422.xlsx")
ext_200bp <- read_xlsx("ws286_ext_Gene3pUtrTrans_200bp_20230422.xlsx")
ext_250bp <- read_xlsx("ws286_ext_Gene3pUtrTrans_250bp_20230422.xlsx")
ext_300bp <- read_xlsx("ws286_ext_Gene3pUtrTrans_300bp_20230422.xlsx")
ext_400bp <- read_xlsx("ws286_ext_Gene3pUtrTrans_400bp_20230422.xlsx")
ext_500bp <- read_xlsx("ws286_ext_Gene3pUtrTrans_500bp_20230422.xlsx")

raw %<>% filter(feature == "gene" & gene_id %in% hh$gene_id) %>% arrange(seqname, gene_id) %>% select(c("gene_id", "start", "end"))
ext_50bp %<>% filter(feature == "gene" & gene_id %in% hh$gene_id) %>% arrange(seqname, gene_id) %>% select(c("gene_id", "start", "end"))
ext_100bp %<>% filter(feature == "gene" & gene_id %in% hh$gene_id) %>% arrange(seqname, gene_id) %>% select(c("gene_id", "start", "end"))
ext_150bp %<>% filter(feature == "gene" & gene_id %in% hh$gene_id) %>% arrange(seqname, gene_id) %>% select(c("gene_id", "start", "end"))
ext_200bp %<>% filter(feature == "gene" & gene_id %in% hh$gene_id) %>% arrange(seqname, gene_id) %>% select(c("gene_id", "start", "end"))
ext_250bp %<>% filter(feature == "gene" & gene_id %in% hh$gene_id) %>% arrange(seqname, gene_id) %>% select(c("gene_id", "start", "end"))
ext_300bp %<>% filter(feature == "gene" & gene_id %in% hh$gene_id) %>% arrange(seqname, gene_id) %>% select(c("gene_id", "start", "end"))
ext_400bp %<>% filter(feature == "gene" & gene_id %in% hh$gene_id) %>% arrange(seqname, gene_id) %>% select(c("gene_id", "start", "end"))
ext_500bp %<>% filter(feature == "gene" & gene_id %in% hh$gene_id) %>% arrange(seqname, gene_id) %>% select(c("gene_id", "start", "end"))

hh %<>% unique()

before_discovering_dcc_has_weird_sorting_order <- cbind(raw,
                                                        ext_50bp[, 2:3],
                                                        ext_100bp[, 2:3],
                                                        ext_150bp[, 2:3],
                                                        ext_200bp[, 2:3],
                                                        ext_250bp[, 2:3],
                                                        ext_300bp[, 2:3],
                                                        ext_400bp[, 2:3],
                                                        ext_500bp[, 2:3])

colnames(before_discovering_dcc_has_weird_sorting_order)[4:19] <- c("start_50bp",
                                                                    "end_50bp",
                                                                    "start_100bp",
                                                                    "end_100bp",
                                                                    "start_150bp",
                                                                    "end_150bp",
                                                                    "start_200bp",
                                                                    "end_200bp",
                                                                    "start_250bp",
                                                                    "end_250bp",
                                                                    "start_300bp",
                                                                    "end_300bp",
                                                                    "start_400bp",
                                                                    "end_400bp",
                                                                    "start_500bp",
                                                                    "end_500bp")

hh_long <- hh %>% pivot_longer(cols = -c(gene_id, strand, end, end_50bp, end_100bp, end_150bp, end_200bp, end_250bp, end_300bp, end_400bp, end_500bp), names_to = "ext", values_to = "start")
hh_long %<>% pivot_longer(cols = -c(gene_id, strand, start, ext), names_to = "ext_end", values_to = "end")
hh_long$ext %<>%
  str_replace("start_50bp", "50") %>%
  str_replace("start_100bp", "100") %>%
  str_replace("start_150bp", "150") %>%
  str_replace("start_200bp", "200") %>%
  str_replace("start_250bp", "250") %>%
  str_replace("start_300bp", "300") %>%
  str_replace("start_400bp", "400") %>%
  str_replace("start_500bp", "500") %>%
  str_replace("start", "0")

hh_long$ext_end %<>%
  str_replace("end_50bp", "50") %>%
  str_replace("end_100bp", "100") %>%
  str_replace("end_150bp", "150") %>%
  str_replace("end_200bp", "200") %>%
  str_replace("end_250bp", "250") %>%
  str_replace("end_300bp", "300") %>%
  str_replace("end_400bp", "400") %>%
  str_replace("end_500bp", "500") %>%
  str_replace("end", "0")

hh_long <- hh_long[hh_long$ext == hh_long$ext_end, ]
hh_long %<>% select(-c("ext_end"))

before_disc_long <- before_discovering_dcc_has_weird_sorting_order %>% pivot_longer(cols = -c(gene_id, end, end_50bp, end_100bp, end_150bp, end_200bp, end_250bp, end_300bp, end_400bp, end_500bp), names_to = "ext", values_to = "start")
before_disc_long %<>% pivot_longer(cols = -c(gene_id, start, ext), names_to = "ext_end", values_to = "end")
before_disc_long$ext %<>%
  str_replace("start_50bp", "50") %>%
  str_replace("start_100bp", "100") %>%
  str_replace("start_150bp", "150") %>%
  str_replace("start_200bp", "200") %>%
  str_replace("start_250bp", "250") %>%
  str_replace("start_300bp", "300") %>%
  str_replace("start_400bp", "400") %>%
  str_replace("start_500bp", "500") %>%
  str_replace("start", "0")

before_disc_long$ext_end %<>%
  str_replace("end_50bp", "50") %>%
  str_replace("end_100bp", "100") %>%
  str_replace("end_150bp", "150") %>%
  str_replace("end_200bp", "200") %>%
  str_replace("end_250bp", "250") %>%
  str_replace("end_300bp", "300") %>%
  str_replace("end_400bp", "400") %>%
  str_replace("end_500bp", "500") %>%
  str_replace("end", "0")

before_disc_long <- before_disc_long[before_disc_long$ext == before_disc_long$ext_end, ]
before_disc_long %<>% select(-c("ext_end"))

colnames(before_disc_long)[3:4] <- c("start_before_disc", "end_before_disc")

before_disc_long$ext %<>% as.numeric()
before_disc_long %<>% arrange(gene_id, ext)

hh_long$ext %<>% as.numeric()
hh_long %<>% arrange(gene_id, ext)

comp <- cbind(hh_long, before_disc_long[, 3:4])
comp$equal <- NA

for (i in 1:nrow(comp)) {
  if(comp$strand[i] == "+") {
    comp$equal[i] <- (comp$end[i] == comp$end_before_disc[i])
  } else {
    comp$equal[i] <- (comp$start[i] == comp$start_before_disc[i])
  }
}

# unbelievable
# before realizing DCC conda R uses the opposite rule of sorting
# although I had to manually change the edits for three genes
# everything else in the extensions was correct
which(!(comp$equal))
