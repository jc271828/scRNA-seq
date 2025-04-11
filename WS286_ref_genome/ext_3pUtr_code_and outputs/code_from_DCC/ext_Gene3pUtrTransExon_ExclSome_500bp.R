library(tidyverse)
library(magrittr)
library(readxl)
library(writexl)

setwd("/work/jc750/ref_genome/WS286/extend_3pUtr")

# load data
spillover <- read_xlsx("./ws286_find_spillovers_ShouldYouExclSome_400bp.xlsx", col_names = T) %>% filter((StopCoordSpillOver == TRUE) & (should_exclude == FALSE)) %>% .$gene_id

# define function to ext gene stop coordinates per chr
# df has to be pre-sorted in ascending order of start position FOR EACH CHROMOSOME!
ExtGenePerChr <- function(chr_df, chr_max, ext) {
  chr_df$start_new <- NA
  chr_df$end_new <- NA
  n <- nrow(chr_df)
  
  # extend all rows except for 1st and last rows
  chr_df[n+1, "start"] <- chr_max + 1
  chr_df[n+1, "end"] <- chr_max + 1
  chr_df[n+1, "strand"] <- chr_df[n, "strand"]
  
  chr_df %<>% rbind(NA, .)
  chr_df[1, "start"] <- 0
  chr_df[1, "end"] <- 0
  chr_df[1, "strand"] <- chr_df[2, "strand"]
  
  for (i in 2:(n+1)) {
    if (chr_df$gene_id[i] %in% spillover) {
      chr_df$start_new[i] <- chr_df$start[i]
      chr_df$end_new[i] <- chr_df$end[i]
    } else {
      if (chr_df$strand[i] == "+") {
        chr_df$start_new[i] <- chr_df$start[i]
        next_gene <- chr_df %>% filter(start > chr_df$end[i]) %>% arrange(start, desc(strand)) %>% .[1, ]
        if (next_gene$strand == "+") {
          if (chr_df$end[i] + ext < next_gene$start) {
            chr_df$end_new[i] <- chr_df$end[i] + ext
          } else {
            chr_df$end_new[i] <- next_gene$start - 1
          }
        } else {
          if (chr_df$end[i] + (2 * ext) < next_gene$start) {
            chr_df$end_new[i] <- chr_df$end[i] + ext
          } else {
            dis <- next_gene$start - chr_df$end[i]
            chr_df$end_new[i] <- chr_df$end[i] + floor(dis/2)
          }
        }
      } else {
        chr_df$end_new[i] <- chr_df$end[i]
        previous_gene <- chr_df %>% filter(end < chr_df$start[i]) %>% arrange(desc(end), strand) %>% .[1, ]
        if (previous_gene$strand == "-") {
          if (chr_df$start[i] - ext  > previous_gene$end) {
            chr_df$start_new[i] <- chr_df$start[i] - ext
          } else {
            chr_df$start_new[i] <- previous_gene$end + 1
          }
        } else {
          if (chr_df$start[i] - (2 * ext) > previous_gene$end) {
            chr_df$start_new[i] <- chr_df$start[i] - ext
          } else {
            dis <- chr_df$start[i] - previous_gene$end
            chr_df$start_new[i] <- chr_df$start[i] - floor(dis/2)
          }
        }
      }
    }
  }
  chr_df <- chr_df[2:(n+1), c(1:5, 21:22, 6:20)]
  return(chr_df)
}

# define function to ext gene stop coordinates for an annotation df
ExtGene <- function(df_to_edit, ext) {
  # extract gene features by chromosome
  edit_I <- df_to_edit %>% filter(feature == "gene" & seqname == "I") %>% arrange(start)
  edit_II <- df_to_edit %>% filter(feature == "gene" & seqname == "II") %>% arrange(start)
  edit_III <- df_to_edit %>% filter(feature == "gene" & seqname == "III") %>% arrange(start)
  edit_IV <- df_to_edit %>% filter(feature == "gene" & seqname == "IV") %>% arrange(start)
  edit_V <- df_to_edit %>% filter(feature == "gene" & seqname == "V") %>% arrange(start)
  edit_X <- df_to_edit %>% filter(feature == "gene" & seqname == "X") %>% arrange(start)
  edit_MtDNA <- df_to_edit %>% filter(feature == "gene" & seqname == "MtDNA") %>% arrange(start)
  
  # chr_max
  # I 15072434
  # II 15279421
  # III 13783801
  # IV 17493829
  # V 20924180
  # X 17718942
  # MtDNA 13794 (circular)
  
  edit_I %<>% ExtGenePerChr(15072434, ext)
  edit_II %<>% ExtGenePerChr(15279421, ext)
  edit_III %<>% ExtGenePerChr(13783801, ext)
  edit_IV %<>% ExtGenePerChr(17493829, ext)
  edit_V %<>% ExtGenePerChr(20924180, ext)
  edit_X %<>% ExtGenePerChr(17718942, ext)
  edit_MtDNA %<>% ExtGenePerChr(13794, ext)
  
  df_gene_edited <- rbind(edit_I, edit_II, edit_III, edit_IV, edit_V, edit_X, edit_MtDNA)
  
  df_non_gene <- df_to_edit
  df_non_gene$start_new <- NA
  df_non_gene$end_new <- NA
  df_non_gene <- df_non_gene[, c(1:5, 21:22, 6:20)]
  
  temp <- filter(df_non_gene, feature != "gene") %>% rbind(df_gene_edited) %>% arrange(seqname, gene_id, feature)
  
  return(temp)
}

# define function to ext 3pUtr, Transcript, and Exon
Ext3pUtrTransExon <- function(df) {
  df %<>% arrange(seqname, gene_id, start)
  df$row_num <- rownames(df)
  all_gene_id <- df$gene_id %>% unique()
  
  for (i in all_gene_id) {
    temp <- filter(df, gene_id == i)
    which_strand <- temp$strand %>% unique()
    if (which_strand == "+") {
      old_gene_stop <- filter(temp, feature == "gene") %>% .$end
      new_gene_stop <- filter(temp, feature == "gene") %>% .$end_new
      utr_row <- temp %>% filter(feature == "three_prime_utr" & end == old_gene_stop) %>% .$row_num
      trans_row <- temp %>% filter(feature == "transcript" & end == old_gene_stop) %>% .$row_num
      exon_row <- temp %>% filter(feature == "exon" & end == old_gene_stop) %>% .$row_num
      df[c(utr_row, trans_row, exon_row), "end_new"] <- new_gene_stop
    } else {
      old_gene_stop <- filter(temp, feature == "gene") %>% .$start
      new_gene_stop <- filter(temp, feature == "gene") %>% .$start_new
      utr_row <- temp %>% filter(feature == "three_prime_utr" & start == old_gene_stop) %>% .$row_num
      trans_row <- temp %>% filter(feature == "transcript" & start == old_gene_stop) %>% .$row_num
      exon_row <- temp %>% filter(feature == "exon" & start == old_gene_stop) %>% .$row_num
      df[c(utr_row, trans_row, exon_row), "start_new"] <- new_gene_stop
    }
  }
  
  df[is.na(df$start_new), "start_new"] <- df[is.na(df$start_new), "start"]
  df[is.na(df$end_new), "end_new"] <- df[is.na(df$end_new), "end"]
  
  df %<>% select(-c(start, end, row_num))
  colnames(df)[4:5] <- c("start", "end")
  
  return(df)
}

# do 400bp -> 500bp gene extension followed by 3pUtr, Transcript, Exon extension
read_xlsx("./ws286_ext_Gene3pUtrTrans_400bp_20230422.xlsx", col_names = T) %>%
  arrange(seqname, start) %>%
  ExtGene(100) %>%
  Ext3pUtrTransExon() %>%
  write_xlsx("ws286_ext_Gene3pUtrTrans_500bp_20230422.xlsx", col_names = T)
