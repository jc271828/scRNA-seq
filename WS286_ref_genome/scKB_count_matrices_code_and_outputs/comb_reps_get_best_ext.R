library(tidyverse)
library(magrittr)
library(cowplot)
library(writexl)

# define function to load and combine data from the same reps
ext_mtx_and_comb <- function(x) {
  pat <- paste("_", x, sep = "")
  
  matrices <- list.files() %>% .[str_detect(., pat)]
  for (i in 1:length(matrices)){
    assign(str_remove_all(matrices[i], "\\.csv"), read_csv(matrices[i], col_names = T))
  }
  
  df <- ls()[str_detect(ls(), pat)]
  dflist <- list()
  for (i in 1:length(df)) {
    dflist[[i]] <- eval(as.symbol(df[i]))
  }
  names(dflist) <- df
  
  for (i in 1:length(df)) {
    dflist[[i]]$gene_id %<>% str_remove_all("\\.")
    colnames(dflist[[i]])[2:4] %<>% paste(df[i], sep = "_")
  }
  
  list2env(dflist, envir = .GlobalEnv)
  
  temp <- bind_cols(dflist)
  temp <- cbind(temp[, 1], temp[, str_detect(colnames(temp), "gene_id", negate = T)])
  colnames(temp)[1] <- "gene_id"
  
  return(temp)
}

sckb_ws286_raw <- ext_mtx_and_comb("raw")
sckb_ws286_50bp <- ext_mtx_and_comb("50bp")
sckb_ws286_100bp <- ext_mtx_and_comb("100bp")
sckb_ws286_150bp <- ext_mtx_and_comb("150bp")
sckb_ws286_200bp <- ext_mtx_and_comb("200bp")
sckb_ws286_250bp <- ext_mtx_and_comb("250bp")
sckb_ws286_300bp <- ext_mtx_and_comb("300bp")
sckb_ws286_400bp <- ext_mtx_and_comb("400bp")
sckb_ws286_500bp <- ext_mtx_and_comb("500bp")

rm(list = ls()[str_detect(ls(), "_sckb_mtx|_scKB")])

# how variable different reps are regarding the same condition (fed/starv) and same extension (only considering sp + unsp counts)
total_counts_per_rep <- function(df) {
  temp_fed <- df[, str_detect(colnames(df), "total_fed")]
  temp_starv <- df[, str_detect(colnames(df), "total_starv")]
  df <- cbind(df[, 1], temp_fed, temp_starv)
  df$total_count_fed_each_lib <- apply(temp_fed, MARGIN = 1, FUN = function(x) x %>% as.numeric() %>% paste(collapse = ","))
  df$total_count_starv_each_lib <- apply(temp_starv, MARGIN = 1, FUN = function(x) x %>% as.numeric() %>% paste(collapse = ","))
  colnames(df)[1] <- "gene_id"
  return(df)
}

# whether an extension gains or loses reads varies by rep
# so it seems like a bad idea to use only one rep to determine the best extension
count_per_lib_raw <- total_counts_per_rep(sckb_ws286_raw)
count_per_lib_50bp <- total_counts_per_rep(sckb_ws286_50bp)
count_per_lib_100bp <- total_counts_per_rep(sckb_ws286_100bp)
count_per_lib_150bp <- total_counts_per_rep(sckb_ws286_150bp)
count_per_lib_200bp <- total_counts_per_rep(sckb_ws286_200bp)
count_per_lib_250bp <- total_counts_per_rep(sckb_ws286_250bp)
count_per_lib_300bp <- total_counts_per_rep(sckb_ws286_300bp)
count_per_lib_400bp <- total_counts_per_rep(sckb_ws286_400bp)
count_per_lib_500bp <- total_counts_per_rep(sckb_ws286_500bp)

# use BY-REP sp + unsp reads (starv + fed) to determine the best extension per rep
# and only change annotation for genes whose "best" extension is agreed by most reps
total_count <- cbind(count_per_lib_raw[, 1:18],
                     count_per_lib_50bp[, 2:18],
                     count_per_lib_100bp[, 2:18],
                     count_per_lib_150bp[, 2:18],
                     count_per_lib_200bp[, 2:18],
                     count_per_lib_250bp[, 2:18],
                     count_per_lib_300bp[, 2:18],
                     count_per_lib_400bp[, 2:18],
                     count_per_lib_500bp[, 2:18])

total_count <- total_count[, str_detect(colnames(total_count), "4-1", negate = T)]
total_count_rep0 <- total_count[, str_detect(colnames(total_count), "gene_id|fed0|starv0")]
total_count_rep1 <- total_count[, str_detect(colnames(total_count), "gene_id|fed1|starv1")]
total_count_rep2 <- total_count[, str_detect(colnames(total_count), "gene_id|fed2|starv2")]
total_count_rep3 <- total_count[, str_detect(colnames(total_count), "gene_id|fed3|starv3")]
total_count_rep4.2 <- total_count[, str_detect(colnames(total_count), "gene_id|fed4-2|starv4-2")]
total_count_rep5 <- total_count[, str_detect(colnames(total_count), "gene_id|fed5|starv5")]
total_count_rep6.1 <- total_count[, str_detect(colnames(total_count), "gene_id|fed6-1|starv6-1")]
total_count_rep6.2 <- total_count[, str_detect(colnames(total_count), "gene_id|fed6-2|starv6-2")]

# plot CPM distribution
plot_cpm <- function(df, col_num) {
  df[[col_num]] <- df[[col_num]]/sum(df[[col_num]])*(10^6)

  names <- str_split(colnames(df)[col_num], "_") %>% unlist()
  cond <- names %>% .[str_detect(., "fed|starv")]
  ext <- names %>% .[str_detect(., "raw|bp")]

  p <-
    ggplot(df, aes(df[[col_num]])) +
    stat_ecdf() +
    theme_classic() +
    theme(aspect.ratio = 1) +
    coord_cartesian(xlim = c(0, 20)) +
    scale_x_continuous(expand = expansion(mult = c(0, 0)),
                       breaks = seq(from = 0, to = 20, by = 5)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0)),
                       breaks = seq(from = 0, to = 1, by = 0.1)) +
    labs(title = paste(cond, ext, sep = "_"),
         x = "CPM")

  return(p)
}

plot_cpm_per_rep <- function(df) {
  plot_grid(plot_cpm(df, 2),
            plot_cpm(df, 3),
            plot_cpm(df, 4),
            plot_cpm(df, 5),
            plot_cpm(df, 6),
            plot_cpm(df, 7),
            plot_cpm(df, 8),
            plot_cpm(df, 9),
            plot_cpm(df, 10),
            plot_cpm(df, 11),
            plot_cpm(df, 12),
            plot_cpm(df, 13),
            plot_cpm(df, 14),
            plot_cpm(df, 15),
            plot_cpm(df, 16),
            plot_cpm(df, 17),
            plot_cpm(df, 18),
            plot_cpm(df, 19),
            nrow = 5,
            ncol = 4)
}

# save CPM distributions to a PDF
cairo_pdf(onefile = T, filename = "CPM_distr_PerRepPerCondPerExt.pdf", width = 10, height = 10)
plot_cpm_per_rep(total_count_rep0)
plot_cpm_per_rep(total_count_rep1)
plot_cpm_per_rep(total_count_rep2)
plot_cpm_per_rep(total_count_rep3)
plot_cpm_per_rep(total_count_rep4.2)
plot_cpm_per_rep(total_count_rep5)
plot_cpm_per_rep(total_count_rep6.1)
plot_cpm_per_rep(total_count_rep6.2)
dev.off()

# use CPM of fed + starv instead of raw read counts of fed + starv to determine best extension
# otherwise, if fed and starv (same rep) sequencing depth differs significantly,
# the deeper library is going to dominate the picking
convert_count_to_cpm <- function(df) {
  for (i in 2: ncol(df)) {
    df[[i]] <- df[[i]]/sum(df[[i]])*(10^6)
  }
  return(df)
}
total_count_rep0 %<>% convert_count_to_cpm()
total_count_rep1 %<>% convert_count_to_cpm()
total_count_rep2 %<>% convert_count_to_cpm()
total_count_rep3 %<>% convert_count_to_cpm()
total_count_rep4.2 %<>% convert_count_to_cpm()
total_count_rep5 %<>% convert_count_to_cpm()
total_count_rep6.1 %<>% convert_count_to_cpm()
total_count_rep6.2 %<>% convert_count_to_cpm()

# get fed + starv CPM per gene per extension for each rep
# underscore is SUPER important in detection pattern. Without it, 150bp, 250bp will also be detected as 50bp!!
cal_fed_plus_starv <- function(df) {
  df$total_raw <- df[, str_detect(colnames(df), "_raw")] %>% rowSums()
  df$total_50bp <- df[, str_detect(colnames(df), "_50bp")] %>% rowSums()
  df$total_100bp <- df[, str_detect(colnames(df), "_100bp")] %>% rowSums()
  df$total_150bp <- df[, str_detect(colnames(df), "_150bp")] %>% rowSums()
  df$total_200bp <- df[, str_detect(colnames(df), "_200bp")] %>% rowSums()
  df$total_250bp <- df[, str_detect(colnames(df), "_250bp")] %>% rowSums()
  df$total_300bp <- df[, str_detect(colnames(df), "_300bp")] %>% rowSums()
  df$total_400bp <- df[, str_detect(colnames(df), "_400bp")] %>% rowSums()
  df$total_500bp <- df[, str_detect(colnames(df), "_500bp")] %>% rowSums()
  return(df)
}

total_count_rep0 %<>% cal_fed_plus_starv()
total_count_rep1 %<>% cal_fed_plus_starv()
total_count_rep2 %<>% cal_fed_plus_starv()
total_count_rep3 %<>% cal_fed_plus_starv()
total_count_rep4.2 %<>% cal_fed_plus_starv()
total_count_rep5 %<>% cal_fed_plus_starv()
total_count_rep6.1 %<>% cal_fed_plus_starv()
total_count_rep6.2 %<>% cal_fed_plus_starv()

# calculate the largest number of reads gained by extension and which extension that gain corresponds to
cal_max_gain_and_percent_gain <- function(df) {
  df$max_CPM <- apply(df[, 20:28], MARGIN = 1, max)
  df$max_gain <- df$max_CPM - df$total_raw
  df$max_gain_ratio_baseline <- NA
  df$percent_gain_50bp <- NA
  df$percent_gain_100bp <- NA
  df$percent_gain_150bp <- NA
  df$percent_gain_200bp <- NA
  df$percent_gain_250bp <- NA
  df$percent_gain_300bp <- NA
  df$percent_gain_400bp <- NA
  df$percent_gain_500bp <- NA
  
  perc_col_num <- which(str_detect(colnames(df), "percent"))
  
  for (i in 1:nrow(df)) {
    if (df[i, "max_gain"] == 0) {
      df[i, perc_col_num] <- 0
    } else {
      df[i, "percent_gain_50bp"] <- (df[i, "total_50bp"] - df[i, "total_raw"])/df[i, "max_gain"]
      df[i, "percent_gain_100bp"] <- (df[i, "total_100bp"] - df[i, "total_raw"])/df[i, "max_gain"]
      df[i, "percent_gain_150bp"] <- (df[i, "total_150bp"] - df[i, "total_raw"])/df[i, "max_gain"]
      df[i, "percent_gain_200bp"] <- (df[i, "total_200bp"] - df[i, "total_raw"])/df[i, "max_gain"]
      df[i, "percent_gain_250bp"] <- (df[i, "total_250bp"] - df[i, "total_raw"])/df[i, "max_gain"]
      df[i, "percent_gain_300bp"] <- (df[i, "total_300bp"] - df[i, "total_raw"])/df[i, "max_gain"]
      df[i, "percent_gain_400bp"] <- (df[i, "total_400bp"] - df[i, "total_raw"])/df[i, "max_gain"]
      df[i, "percent_gain_500bp"] <- (df[i, "total_500bp"] - df[i, "total_raw"])/df[i, "max_gain"]
    }
    
    if (df[i, "total_raw"] == 0) {
      df[i, "max_gain_ratio_baseline"] <- NA
    } else {
      df[i, "max_gain_ratio_baseline"] <- df[i, "max_gain"]/df[i, "total_raw"]
    }
  }
  
  return(df)
}

total_count_rep0 %<>% cal_max_gain_and_percent_gain()
total_count_rep1 %<>% cal_max_gain_and_percent_gain()
total_count_rep2 %<>% cal_max_gain_and_percent_gain()
total_count_rep3 %<>% cal_max_gain_and_percent_gain()
total_count_rep4.2 %<>% cal_max_gain_and_percent_gain()
total_count_rep5 %<>% cal_max_gain_and_percent_gain()
total_count_rep6.1 %<>% cal_max_gain_and_percent_gain()
total_count_rep6.2 %<>% cal_max_gain_and_percent_gain()

# get an idea of how max_CPM distributes (need to pick a cutoff: max_CPM >= 2 or >= 1 or else?)
# pick max_CPM >= 1
plot_max_CPM <- function(df) {
  p <- 
    ggplot(df, aes(max_CPM)) +
    stat_ecdf() +
    theme_classic() +
    theme(aspect.ratio = 1) +
    coord_cartesian(xlim = c(0, 20)) +
    scale_x_continuous(expand = expansion(mult = c(0, 0)),
                       breaks = seq(from = 0, to = 20, by = 1)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0)),
                       breaks = seq(from = 0, to = 1, by = 0.1)) +
    labs(title = deparse(substitute(df)))
  
  return(p)
}

plot_grid(plot_max_CPM(total_count_rep0),
          plot_max_CPM(total_count_rep1),
          plot_max_CPM(total_count_rep2),
          plot_max_CPM(total_count_rep3),
          plot_max_CPM(total_count_rep4.2),
          plot_max_CPM(total_count_rep5),
          plot_max_CPM(total_count_rep6.1),
          plot_max_CPM(total_count_rep6.2),
          nrow = 3,
          ncol = 3)

# mark max_CPM <1 for removal
total_count_rep0$max_CPM_LessThan1 <- (total_count_rep0$max_CPM < 1)
total_count_rep1$max_CPM_LessThan1 <- (total_count_rep1$max_CPM < 1)
total_count_rep2$max_CPM_LessThan1 <- (total_count_rep2$max_CPM < 1)
total_count_rep3$max_CPM_LessThan1 <- (total_count_rep3$max_CPM < 1)
total_count_rep4.2$max_CPM_LessThan1 <- (total_count_rep4.2$max_CPM < 1)
total_count_rep5$max_CPM_LessThan1 <- (total_count_rep5$max_CPM < 1)
total_count_rep6.1$max_CPM_LessThan1 <- (total_count_rep6.1$max_CPM < 1)
total_count_rep6.2$max_CPM_LessThan1 <- (total_count_rep6.2$max_CPM < 1)

# how many genes have past the CPM>=1 filter in each rep to have a chance of being modified (have to pass other filters too though)?
total_count_rep0 %>% filter(max_CPM_LessThan1 == FALSE) %>% dim() # 17433
total_count_rep1 %>% filter(max_CPM_LessThan1 == FALSE) %>% dim() # 14799
total_count_rep2 %>% filter(max_CPM_LessThan1 == FALSE) %>% dim() # 15509
total_count_rep3 %>% filter(max_CPM_LessThan1 == FALSE) %>% dim() # 14424
total_count_rep4.2 %>% filter(max_CPM_LessThan1 == FALSE) %>% dim() # 14752
total_count_rep5 %>% filter(max_CPM_LessThan1 == FALSE) %>% dim() # 15569
total_count_rep6.1 %>% filter(max_CPM_LessThan1 == FALSE) %>% dim() # 14843
total_count_rep6.2 %>% filter(max_CPM_LessThan1 == FALSE) %>% dim() # 14678

# get an idea of how max_gain_ratio_baseline distributes (need to pick a cutoff)
# pick max_gain_ratio_baseline >= 0.05
plot_max_gain <- function(df) {
  p <- 
    ggplot(df, aes(max_gain_ratio_baseline)) +
    stat_ecdf() +
    theme_classic() +
    theme(aspect.ratio = 1) +
    coord_cartesian(xlim = c(0, 0.25)) +
    scale_x_continuous(expand = expansion(mult = c(0, 0))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0)),
                       breaks = seq(from = 0, to = 1, by = 0.1)) +
    labs(title = deparse(substitute(df)))
  
  return(p)
}

plot_max_gain_CPM_filtered <- function(df) {
  p <- 
    ggplot(filter(df, max_CPM_LessThan1 == FALSE), aes(max_gain_ratio_baseline)) +
    stat_ecdf() +
    theme_classic() +
    theme(aspect.ratio = 1) +
    coord_cartesian(xlim = c(0, 0.25)) +
    scale_x_continuous(expand = expansion(mult = c(0, 0))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0)),
                       breaks = seq(from = 0, to = 1, by = 0.1)) +
    labs(title = deparse(substitute(df)))
  
  return(p)
}

plot_grid(plot_max_gain_CPM_filtered(total_count_rep0),
          plot_max_gain_CPM_filtered(total_count_rep1),
          plot_max_gain_CPM_filtered(total_count_rep2),
          plot_max_gain_CPM_filtered(total_count_rep3),
          plot_max_gain_CPM_filtered(total_count_rep4.2),
          plot_max_gain_CPM_filtered(total_count_rep5),
          plot_max_gain_CPM_filtered(total_count_rep6.1),
          plot_max_gain_CPM_filtered(total_count_rep6.2),
          nrow = 3,
          ncol = 3)

# quantile(ecdf(total_count_rep0$max_gain_ratio_baseline), 0.95)
# quantile(ecdf(total_count_rep1$max_gain_ratio_baseline), 0.95)
# quantile(ecdf(total_count_rep2$max_gain_ratio_baseline), 0.95)
# quantile(ecdf(total_count_rep3$max_gain_ratio_baseline), 0.95)
# quantile(ecdf(total_count_rep4.2$max_gain_ratio_baseline), 0.95)
# quantile(ecdf(total_count_rep5$max_gain_ratio_baseline), 0.95)
# quantile(ecdf(total_count_rep6.1$max_gain_ratio_baseline), 0.95)
# quantile(ecdf(total_count_rep6.2$max_gain_ratio_baseline), 0.95)

# pick the "best" extension for each gene in each rep
# "best": shortest extension that gains at least 90% of max gain
pick_shortest_ext_with_90p_max_gain <- function(df) {
  df$best_ext <- NA
  
  temp <- filter(df, max_CPM_LessThan1 == FALSE & (max_gain_ratio_baseline >= 0.05 | total_raw == 0))
  perc_col_num <- which(str_detect(colnames(temp), "percent"))
  
  for (i in 1:nrow(temp)) {
    AtLeast90pMaxGain_col_num <- perc_col_num[temp[i, perc_col_num] >= 0.9]
    temp[i, "best_ext"] <- AtLeast90pMaxGain_col_num %>% min() %>% colnames(df)[.] %>% str_split("_|bp") %>% unlist() %>% .[3] %>% as.numeric()
  }
  
  df %<>% filter(!(max_CPM_LessThan1 == FALSE & (max_gain_ratio_baseline >= 0.05 | total_raw == 0))) %>% rbind(temp) %>% arrange(gene_id)
  df$best_ext[is.na(df$best_ext)] <- 0
  
  return(df)
}

total_count_rep0 %<>% pick_shortest_ext_with_90p_max_gain()
total_count_rep1 %<>% pick_shortest_ext_with_90p_max_gain()
total_count_rep2 %<>% pick_shortest_ext_with_90p_max_gain()
total_count_rep3 %<>% pick_shortest_ext_with_90p_max_gain()
total_count_rep4.2 %<>% pick_shortest_ext_with_90p_max_gain()
total_count_rep5 %<>% pick_shortest_ext_with_90p_max_gain()
total_count_rep6.1 %<>% pick_shortest_ext_with_90p_max_gain()
total_count_rep6.2 %<>% pick_shortest_ext_with_90p_max_gain()

colnames(total_count_rep0)[41] %<>% paste("rep0", sep = "_")
colnames(total_count_rep1)[41] %<>% paste("rep1", sep = "_")
colnames(total_count_rep2)[41] %<>% paste("rep2", sep = "_")
colnames(total_count_rep3)[41] %<>% paste("rep3", sep = "_")
colnames(total_count_rep4.2)[41] %<>% paste("rep4.2", sep = "_")
colnames(total_count_rep5)[41] %<>% paste("rep5", sep = "_")
colnames(total_count_rep6.1)[41] %<>% paste("rep6.1", sep = "_")
colnames(total_count_rep6.2)[41] %<>% paste("rep6.2", sep = "_")

# combine all reps and see how well reps agree on the best extension length
best_ext_all_reps <- cbind(total_count_rep0[, c(1, 41)],
                           total_count_rep1$best_ext_rep1,
                           total_count_rep2$best_ext_rep2,
                           total_count_rep3$best_ext_rep3,
                           total_count_rep4.2$best_ext_rep4.2,
                           total_count_rep5$best_ext_rep5,
                           total_count_rep6.1$best_ext_rep6.1,
                           total_count_rep6.2$best_ext_rep6.2)

colnames(best_ext_all_reps)[3:9] <- c("best_ext_rep1",
                                      "best_ext_rep2",
                                      "best_ext_rep3",
                                      "best_ext_rep4.2",
                                      "best_ext_rep5",
                                      "best_ext_rep6.1",
                                      "best_ext_rep6.2")

best_ext_all_reps$best_ext <- NA
best_ext_all_reps$num_reps_agreed <- NA

# at least 2 reps agree on an extension -> can pick it
for (i in 1:nrow(best_ext_all_reps)) {
  exts <- best_ext_all_reps[i, 2:9] %>% as.numeric() %>% table() %>% as.data.frame()
  colnames(exts)[1] <- "Ext"
  exts$Ext %<>% as.character() %>% as.numeric() # somehow have to convert to character and then number
  exts %<>% arrange(desc(Freq), desc(Ext))
  best_ext_all_reps$best_ext[i] <- exts[1, 1] %>% as.numeric()
  best_ext_all_reps$num_reps_agreed[i] <- exts[1, 2] %>% as.numeric()
}

best_ext <- best_ext_all_reps %>% filter(num_reps_agreed >= 2 & best_ext != 0) # 1536 genes modified
write_xlsx(best_ext, "BestExt2Reps0.05GainRatio_20230422.xlsx", col_names = T)

# # plot with random gene_id
# cairo_pdf(onefile = T, file = "RandomGenes_AvgReadCount_Per3pUtrExt_CombiningFedAndStarv.pdf", width = 10, height = 10, family = "sans")
# 
# ran_gene <- sckb_ws286_wide$gene_id %>% unique() %>% length() %>% runif(9, min=1, max=.) %>% unique(sckb_ws286_wide$gene_id)[.]
# ggplot(data = filter(sckb_ws286_wide, gene_id %in% ran_gene)) +
#   geom_line(mapping = aes(x = ThreePrimeUtr_extension, y = log2(avg_read_count+1), group = gene_id)) +
#   theme_classic() +
#   theme(aspect.ratio = 1) +
#   facet_wrap(.~gene_id)
# 
# dev.off()
