
# load libraries --------------------------------------------------------

library(tidyverse)
library(magrittr)


# define a function to clean loaded csv files -----------------------------

clean_csv <- function(path) {
  temp <- read_csv(path, col_names = TRUE)
  temp <- temp[, -1]
  colnames(temp)[1] <- "filtered_pgc_id"
  return(temp)
}


# change dir --------------------------------------------------------------

setwd("../../")
parent_dir <- getwd()


# get directories of interest ---------------------------------------------

all_dirs <- list.dirs(path = ".", full.names = TRUE, recursive = FALSE)
matching_dirs <- all_dirs[grep("./202501*", all_dirs)]


# get all GC1171 data frames ----------------------------------------------

GC1171 <- data.frame()

for (i in matching_dirs) {
  setwd(parent_dir)
  setwd(i)
  
  splitted_dir <- str_split(getwd(), "_") %>% unlist()
  
  if (dir.exists("GC1171")) {
    setwd("GC1171")
    
    df_clust_min_3 <- clean_csv("smfish_result_merged_clust_min_3.csv")
    df_clust_min_4 <- clean_csv("smfish_result_merged_clust_min_4.csv")
    df_clust_min_5 <- clean_csv("smfish_result_merged_clust_min_5.csv")
    
    df <- rbind(df_clust_min_3, df_clust_min_4, df_clust_min_5)
    df$rep <- splitted_dir[length(splitted_dir) - 1]
    df$hours_post_bleach <- splitted_dir[length(splitted_dir)]
    
    GC1171 %<>% rbind(df)
  }
}

GC1171 %<>% filter(!(rep == "R2" & hours_post_bleach == "24h"))


# get all LRB684 data frames ----------------------------------------------

LRB684 <- data.frame()

for (i in matching_dirs) {
  setwd(parent_dir)
  setwd(i)
  
  splitted_dir <- str_split(getwd(), "_") %>% unlist()
  
  if (dir.exists("LRB684")) {
    setwd("LRB684")
    
    df_clust_min_3 <- clean_csv("smfish_result_merged_clust_min_3.csv")
    df_clust_min_4 <- clean_csv("smfish_result_merged_clust_min_4.csv")
    df_clust_min_5 <- clean_csv("smfish_result_merged_clust_min_5.csv")
    
    df <- rbind(df_clust_min_3, df_clust_min_4, df_clust_min_5)
    df$rep <- splitted_dir[length(splitted_dir) - 1]
    df$hours_post_bleach <- splitted_dir[length(splitted_dir)]
    
    LRB684 %<>% rbind(df)
  }
}


# clean data --------------------------------------------------------------

GC1171$genotype <- "wildtype"
LRB684$genotype <- "aak-1(del)"

compare <- rbind(GC1171, LRB684)
compare$pgc_has_intronic_signal <- NA

compare$clusters_inside_pgc %>% range()
compare$clusters_inside_pgc %>% table()

compare$pgc_has_intronic_signal[compare$clusters_inside_pgc > 0] <- TRUE
compare$pgc_has_intronic_signal[compare$clusters_inside_pgc == 0] <- FALSE

for_stats <-
  compare %>%
  group_by(genotype, hours_post_bleach, min_spots_per_cluster, pgc_has_intronic_signal) %>%
  summarise(number_of_pgc = n(), .groups = "drop") %>% # .groups = "drop" does the magic. It removes all grouping structures after summarisation. Without it, the number of rows will be wrong
  pivot_wider(names_from = "pgc_has_intronic_signal", values_from = "number_of_pgc") %>%
  arrange(hours_post_bleach, min_spots_per_cluster)

colnames(for_stats)[4:5] <- c("no_ts_in_pgc", "has_ts_in_pgc")
for_stats$total_pgc <- for_stats$no_ts_in_pgc + for_stats$has_ts_in_pgc

# per_rep <-
#   compare %>%
#   group_by(genotype, hours_post_bleach, min_spots_per_cluster, pgc_has_intronic_signal, rep) %>%
#   summarise(number_of_pgc = n(), .groups = "drop") %>% # .groups = "drop" does the magic. It removes all grouping structures after summarisation. Without it, the number of rows will be wrong
#   pivot_wider(names_from = "pgc_has_intronic_signal", values_from = "number_of_pgc") %>%
#   arrange(hours_post_bleach, min_spots_per_cluster, rep)
# 
# colnames(per_rep)[5:6] <- c("no_ts_in_pgc", "has_ts_in_pgc")
# per_rep$total_pgc <- per_rep$no_ts_in_pgc + per_rep$has_ts_in_pgc
# 
# per_rep %<>% filter(hours_post_bleach == "24h")


# do stats ----------------------------------------------------------------

for_stats %<>% mutate(percent_with_ts = round(has_ts_in_pgc/total_pgc * 100, 1))
for_stats %<>% filter(min_spots_per_cluster != 3)

# pick min_spots_per_cluster == 4
for_stats %<>% filter(min_spots_per_cluster == 4)
for_stats %<>% select(-c("total_pgc", "percent_with_ts"))
for_stats %<>% as.data.frame()
row.names(for_stats) <- paste(for_stats$genotype, for_stats$hours_post_bleach, sep = "_")

# manually check images and do corrections
for_stats["wildtype_12h", "no_ts_in_pgc"] <- for_stats["wildtype_12h", "no_ts_in_pgc"] + 1 # _2024_12_28__16_41_33_Out TIFFs\clust_min_4
for_stats["wildtype_12h", "no_ts_in_pgc"] <- for_stats["wildtype_12h", "no_ts_in_pgc"] - 1 # _2025_01_11__12_57_47_Out TIFFs\clust_min_4
for_stats["wildtype_12h", "has_ts_in_pgc"] <- for_stats["wildtype_12h", "has_ts_in_pgc"] + 1 # _2025_01_11__12_57_47_Out TIFFs\clust_min_4
for_stats["wildtype_12h", "has_ts_in_pgc"] <- for_stats["wildtype_12h", "has_ts_in_pgc"] + 1 # _2025_01_11__13_38_17_Out TIFFs\clust_min_4
for_stats["wildtype_12h", "no_ts_in_pgc"] <- for_stats["wildtype_12h", "no_ts_in_pgc"] + 1 # _2025_01_15__18_15_55_Out TIFFs\clust_min_4
for_stats["wildtype_12h", "no_ts_in_pgc"] <- for_stats["wildtype_12h", "no_ts_in_pgc"] + 1 # _2025_01_15__18_39_03_Out TIFF\clust_min_4
for_stats["wildtype_12h", "has_ts_in_pgc"] <- for_stats["wildtype_12h", "has_ts_in_pgc"] + 1 # _2025_01_15__20_44_04_Out TIFFs\clust_min_4

for_stats["aak-1(del)_12h", "no_ts_in_pgc"] <- for_stats["aak-1(del)_12h", "no_ts_in_pgc"] + 1 # _2024_12_28__19_24_02_Out TIFFs\clust_min_4
for_stats["aak-1(del)_12h", "no_ts_in_pgc"] <- for_stats["aak-1(del)_12h", "no_ts_in_pgc"] + 1 # _2024_12_28__20_38_27_Out TIFFs\clust_min_4
for_stats["aak-1(del)_12h", "no_ts_in_pgc"] <- for_stats["aak-1(del)_12h", "no_ts_in_pgc"] + 1 # _2024_12_28__20_57_40_Out TIFFs\clust_min_4
for_stats["aak-1(del)_12h", "no_ts_in_pgc"] <- for_stats["aak-1(del)_12h", "no_ts_in_pgc"] + 2 # _2025_01_11__15_01_46_Out TIFFs\clust_min_4
for_stats["aak-1(del)_12h", "no_ts_in_pgc"] <- for_stats["aak-1(del)_12h", "no_ts_in_pgc"] + 2 # _2025_01_15__22_50_00_Out TIFFs\clust_min_4

for_stats["wildtype_24h", "no_ts_in_pgc"] <- for_stats["wildtype_24h", "no_ts_in_pgc"] + 1 # _2025_01_10__18_18_43_Out TIFFs\clust_min_4
for_stats["wildtype_24h", "no_ts_in_pgc"] <- for_stats["wildtype_24h", "no_ts_in_pgc"] + 2 # _2025_01_10__18_57_07_Out TIFFs\clust_min_4
for_stats["wildtype_24h", "no_ts_in_pgc"] <- for_stats["wildtype_24h", "no_ts_in_pgc"] + 1 # _2025_01_10__19_38_35_Out TIFFs\clust_min_4
for_stats["wildtype_24h", "no_ts_in_pgc"] <- for_stats["wildtype_24h", "no_ts_in_pgc"] + 1 # _2025_01_10__19_57_57_Out TIFFs\clust_min_4
for_stats["wildtype_24h", "no_ts_in_pgc"] <- for_stats["wildtype_24h", "no_ts_in_pgc"] + 1 # _2025_01_10__20_36_30_Out TIFFs\clust_min_4
for_stats["wildtype_24h", "has_ts_in_pgc"] <- for_stats["wildtype_24h", "has_ts_in_pgc"] + 1 # _2025_01_15__23_13_25_Out TIFFs\clust_min_4
for_stats["wildtype_24h", "no_ts_in_pgc"] <- for_stats["wildtype_24h", "no_ts_in_pgc"] + 1 # _2025_01_15__23_19_23_Out TIFFs\clust_min_4
for_stats["wildtype_24h", "no_ts_in_pgc"] <- for_stats["wildtype_24h", "no_ts_in_pgc"] + 1 # _2025_01_16__00_41_35_Out TIFFs\clust_min_4

for_stats["aak-1(del)_24h", "no_ts_in_pgc"] <- for_stats["aak-1(del)_24h", "no_ts_in_pgc"] + 1 # _2025_01_10__21_28_56_Out TIFFs\clust_min_4
for_stats["aak-1(del)_24h", "no_ts_in_pgc"] <- for_stats["aak-1(del)_24h", "no_ts_in_pgc"] + 1 # _2025_01_10__22_33_54_Out TIFFs\clust_min_4
for_stats["aak-1(del)_24h", "no_ts_in_pgc"] <- for_stats["aak-1(del)_24h", "no_ts_in_pgc"] + 1 # _2025_01_10__23_58_40_Out TIFFs\clust_min_4
for_stats["aak-1(del)_24h", "no_ts_in_pgc"] <- for_stats["aak-1(del)_24h", "no_ts_in_pgc"] + 1 # _2025_01_11__00_25_27_Out TIFFs\clust_min_4
for_stats["aak-1(del)_24h", "no_ts_in_pgc"] <- for_stats["aak-1(del)_24h", "no_ts_in_pgc"] + 1 # _2025_01_11__18_43_24_Out TIFFs\clust_min_4
for_stats["aak-1(del)_24h", "no_ts_in_pgc"] <- for_stats["aak-1(del)_24h", "no_ts_in_pgc"] + 1 # _2025_01_11__19_32_52_Out TIFFs\clust_min_4

# stats for 12h post bleach using at least 4 spots per cluster
contingency_table <- filter(for_stats, hours_post_bleach == "12h" & min_spots_per_cluster == 4)[, c("no_ts_in_pgc", "has_ts_in_pgc")]
contingency_table
fisher.test(contingency_table) # 0.0002129
chisq.test(contingency_table) # 0.0002129

# stats for 24h post bleach using at least 4 spots per cluster
contingency_table <- filter(for_stats, hours_post_bleach == "24h" & min_spots_per_cluster == 4)[, c("no_ts_in_pgc", "has_ts_in_pgc")]
contingency_table
fisher.test(contingency_table) # 0.009783
chisq.test(contingency_table) # 0.01517

for_stats$total_pgc <- for_stats$no_ts_in_pgc + for_stats$has_ts_in_pgc


# states for wild type between 12h and 24h post bleach
contingency_table <- filter(for_stats, genotype == "wildtype" & min_spots_per_cluster == 4)[, c("no_ts_in_pgc", "has_ts_in_pgc")]
contingency_table
fisher.test(contingency_table) # 0.2305
chisq.test(contingency_table) # 0.2404


# another way of doing stats ----------------------------------------------

for_stats_more <- compare %>% select(c("clusters_inside_pgc", "min_spots_per_cluster", "rep", "hours_post_bleach", "genotype")) %>% filter(min_spots_per_cluster == 4)

# manually correct some data
temp <- data.frame(clusters_inside_pgc = c(0, 0, 1, 1, 1),
                   min_spots_per_cluster = rep(4, 5),
                   rep = rep("manual", 5),
                   hours_post_bleach = rep("12h", 5),
                   genotype = rep("wildtype", 5))
for_stats_more %<>% rbind(temp)

temp <- data.frame(clusters_inside_pgc = rep(0, 7),
                   min_spots_per_cluster = rep(4, 7),
                   rep = rep("manual", 7),
                   hours_post_bleach = rep("12h", 7),
                   genotype = rep("aak-1(del)", 7))
for_stats_more %<>% rbind(temp)

temp <- data.frame(clusters_inside_pgc = c(rep(0, 8), 1),
                   min_spots_per_cluster = rep(4, 9),
                   rep = rep("manual", 9),
                   hours_post_bleach = rep("24h", 9),
                   genotype = rep("wildtype", 9))
for_stats_more %<>% rbind(temp)

temp <- data.frame(clusters_inside_pgc = rep(0, 6),
                   min_spots_per_cluster = rep(4, 6),
                   rep = rep("manual", 6),
                   hours_post_bleach = rep("24h", 6),
                   genotype = rep("aak-1(del)", 6))
for_stats_more %<>% rbind(temp)

rm(list = "temp")

for_stats_more$category <- paste(for_stats_more$genotype, for_stats_more$clusters_inside_pgc, sep = "_")
for_stats_more$genotype %<>% as.factor() %>% fct_relevel(c("wildtype", "aak-1(del)"))
for_stats_more$clusters_inside_pgc %<>% as.character()

for_stats_more_to_plot <-
  for_stats_more %>%
  group_by(genotype, hours_post_bleach, clusters_inside_pgc) %>%
  summarise(number_of_pgc = n(), .groups = "drop")

# stats data frames have the same number of entries
all.equal(sum(for_stats$total_pgc), nrow(for_stats_more), sum(for_stats_more_to_plot$number_of_pgc))

# prepare data frame for categorical contingency tables
for_stats_more_to_plot_wider <-
  for_stats_more_to_plot %>%
  pivot_wider(names_from = clusters_inside_pgc, values_from = number_of_pgc)

colnames(for_stats_more_to_plot_wider)[3:5] <- c("pgc_without_ts", "pgc_with_one_ts", "pgc_with_two_ts")
for_stats_more_to_plot_wider$pgc_with_two_ts[is.na(for_stats_more_to_plot_wider$pgc_with_two_ts)] <- 0

# stats for 12h post bleach using at least 4 spots per cluster
contingency_table <- filter(for_stats_more_to_plot_wider, hours_post_bleach == "12h")[, c("pgc_without_ts", "pgc_with_one_ts", "pgc_with_two_ts")]
contingency_table
fisher.test(contingency_table) # 0.0002129
chisq.test(contingency_table) # 0.0006926

# stats for 24h post bleach using at least 4 spots per cluster
contingency_table <- filter(for_stats_more_to_plot_wider, hours_post_bleach == "24h")[, c("pgc_without_ts", "pgc_with_one_ts", "pgc_with_two_ts")]
contingency_table
fisher.test(contingency_table[, 1:2]) # 0.009783
chisq.test(contingency_table[, 1:2]) # 0.01517

# stats for GC1171 using at least 4 spots per cluster
contingency_table <- filter(for_stats_more_to_plot_wider, genotype == "wildtype")[, c("pgc_without_ts", "pgc_with_one_ts", "pgc_with_two_ts")]
contingency_table
fisher.test(contingency_table) # 0.2604
chisq.test(contingency_table) # 0.3206

# stats for LRB684 using at least 4 spots per cluster
contingency_table <- filter(for_stats_more_to_plot_wider, genotype == "aak-1(del)")[, c("pgc_without_ts", "pgc_with_one_ts", "pgc_with_two_ts")]
contingency_table
fisher.test(contingency_table[, 1:2]) # 0.5292
chisq.test(contingency_table[, 1:2]) # 0.6703


# plot --------------------------------------------------------------------

ggplot(for_stats_more_to_plot) +
  geom_col(mapping = aes(x = clusters_inside_pgc, y = number_of_pgc, color = clusters_inside_pgc), fill = NA, linewidth = 1, width = 0.8) +
  scale_color_manual(values = c("black", "red", "darkviolet")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                     breaks = seq(0, 150, by = 25)) +
  theme_classic() +
  theme(aspect.ratio = 1) +
  facet_wrap(hours_post_bleach~genotype)














