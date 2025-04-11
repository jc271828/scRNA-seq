
# load libraries ----------------------------------------------------------

library(tidyverse)
library(magrittr)
library(ggpubr)
library(cowplot)
library(car)


# load data ---------------------------------------------------------------

aak1_12h <- read_csv("../20241120_aak-1_GFP_neg_12h_chr_compaction/processed_auto_3D/chr_compaction_metrics.csv") %>% unique() # picked aak-1/hT2; naSi2 L4s; bleached 5 days later (7 * 4 L4s)
aak1_24h <- read_csv("../20241120_aak-1_GFP_neg_24h_chr_compaction/processed_auto_3D/chr_compaction_metrics.csv") %>% unique() # picked aak-1/hT2; naSi2 L4s; bleached 5 days later (7 * 4 L4s)

WT_12h <- read_csv("../20241120_WT_GFP_neg_12h_chr_compaction/processed_auto_3D/chr_compaction_metrics.csv") %>% unique() # picked +/hT2; naSi2 L4s; bleached 5 days later (7 * 4 L4s)
WT_24h <- read_csv("../20241120_WT_GFP_neg_24h_chr_compaction/processed_auto_3D/chr_compaction_metrics.csv") %>% unique() # picked +/hT2; naSi2 L4s; bleached 5 days later (7 * 4 L4s)

aak1_from_het_moms_12h <- read_csv("../20241127_aak-1_GFP_neg_12h_chr_compaction/processed_auto_3D/chr_compaction_metrics.csv") %>% unique() # picked aak-1/hT2; naSi2 adults in the morning; bleached at late night (1000 adults)
aak1_from_het_moms_24h <- read_csv("../20241127_aak-1_GFP_neg_24h_chr_compaction/processed_auto_3D/chr_compaction_metrics.csv") %>% unique() # picked aak-1/hT2; naSi2 adults in the morning; bleached at late night (1000 adults)

WT_from_het_moms_12h <- read_csv("../20241127_WT_GFP_neg_12h_chr_compaction/processed_auto_3D/chr_compaction_metrics.csv") %>% unique() # picked +/hT2; naSi2 adults in the morning; bleached at late night (1000 adults)
WT_from_het_moms_24h <- read_csv("../20241127_WT_GFP_neg_24h_chr_compaction/processed_auto_3D/chr_compaction_metrics.csv") %>% unique() # picked +/hT2; naSi2 adults in the morning; bleached at late night (1000 adults)


# define a function to clean up df ----------------------------------------

clean_df <- function(df) {
  df <- df[, -1]
  colnames(df)[1:2] <- c("PGC_ID", "slice")
  for (i in 1:nrow(df)) {
    temp <- df$PGC_ID[i] %>% str_split("\\\\") %>% unlist()
    image <- temp[length(temp) - 1]
    image %<>% str_split("_Out") %>% unlist() %>% .[1]
    pgc <- temp[length(temp)]
    df$PGC_ID[i] <- paste(image, pgc, sep = " ")
    
    slice_id <- df$slice[i]
    slice_id %<>% str_split("Out-|Image-|.tif") %>% unlist() %>% .[2]
    df$slice[i] <- slice_id
  }
  df %<>% unique()
  return(df)
}


# define a function to extract best slices --------------------------------

pick_best_slices <- function(df) {
  df %<>% arrange(PGC_ID, slice)
  
  # extract PGCs that have more than 1 slices to pick from
  df_need_to_pick <- df$PGC_ID %>% table() %>% as.data.frame() %>% filter(Freq > 1) %>% .[[1]] %>% as.character() %>% sort()
  df_need_to_pick <- df %>% filter(PGC_ID %in% df_need_to_pick)
  
  # plot correlation of RawIntDen ~ ROI_area
  p <-
    ggplot(data = df_need_to_pick,
           mapping = aes(x = ROI_area, y = RawIntDen)) +
    geom_point(size = 2) +
    geom_smooth(method = "lm") +
    stat_cor(method = "pearson") +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
    theme_classic() +
    theme(aspect.ratio = 1) +
    # facet_wrap(.~PGC_ID) +
    labs(x = "Region of Interest Area",
         y = "Raw Integrated Intensity",
         title = "Each dot is a PGC in a Z-stack slice")
  
  print(p)
  
  # do ROI_area and RawIntDen always agree on the same optimal slice?
  df_need_to_pick_best_by_IntDen <-
    df_need_to_pick %>%
    group_by(PGC_ID) %>%
    arrange(desc(RawIntDen)) %>%
    reframe(best_slice_by_IntDen = slice[1])
  
  df_need_to_pick_best_by_area <-
    df_need_to_pick %>%
    group_by(PGC_ID) %>%
    arrange(desc(ROI_area)) %>%
    reframe(best_slice_by_area = slice[1])
  
  df_need_to_pick_best_slice <- merge(df_need_to_pick_best_by_area, df_need_to_pick_best_by_IntDen, by = "PGC_ID", all = T) %>% arrange(PGC_ID)
  df_need_to_pick_best_slice %<>% mutate(do_they_agree = (best_slice_by_IntDen == best_slice_by_area))
  rm(list = c("df_need_to_pick_best_by_IntDen", "df_need_to_pick_best_by_area"))
  df_need_to_pick_best_slice %<>% arrange(do_they_agree, PGC_ID)
  
  return(df_need_to_pick_best_slice)
}


# generate df of best slices ----------------------------------------------

best_slice_df <- function(df, need_to_pick_df, standard) {
  df_best_slices <- df
  df_best_slices$best_slice <- NA
  df_best_slices$best_slice[df_best_slices$PGC_ID %in% need_to_pick_df$PGC_ID] <- FALSE
  df_best_slices$best_slice[is.na(df_best_slices$best_slice)] <- TRUE
  
  if (standard == "intensity") {
    for(i in which(df_best_slices$best_slice == FALSE)) {
      temp <- need_to_pick_df %>% filter(PGC_ID == df_best_slices$PGC_ID[i]) %>% .$best_slice_by_IntDen
      if(df_best_slices$slice[i] == temp) {
        df_best_slices$best_slice[i] <- TRUE
      }
    }
  } else if (standard == "area") {
    for(i in which(df_best_slices$best_slice == FALSE)) {
      temp <- need_to_pick_df %>% filter(PGC_ID == df_best_slices$PGC_ID[i]) %>% .$best_slice_by_area
      if(df_best_slices$slice[i] == temp) {
        df_best_slices$best_slice[i] <- TRUE
      }
    }
  }
  
  df_best_slices %<>% filter(best_slice == TRUE)
  
  return(df_best_slices)
}


# tidy up data ------------------------------------------------------------

aak1_12h %<>% clean_df()
aak1_24h %<>% clean_df()

WT_12h %<>% clean_df()
WT_24h %<>% clean_df()

aak1_from_het_moms_12h %<>% clean_df()
aak1_from_het_moms_24h %<>% clean_df()

WT_from_het_moms_12h %<>% clean_df()
WT_from_het_moms_24h %<>% clean_df()


# pick optimal slices ----------------

aak1_12h_need_to_pick_best_slice <- aak1_12h %>% pick_best_slices()
aak1_12h_need_to_pick_best_slice$do_they_agree %>% table()

aak1_24h_need_to_pick_best_slice <- aak1_24h %>% pick_best_slices()
aak1_24h_need_to_pick_best_slice$do_they_agree %>% table()

WT_12h_need_to_pick_best_slice <- WT_12h %>% pick_best_slices()
WT_12h_need_to_pick_best_slice$do_they_agree %>% table()

WT_24h_need_to_pick_best_slice <- WT_24h %>% pick_best_slices()
WT_24h_need_to_pick_best_slice$do_they_agree %>% table()

aak1_from_het_moms_12h_need_to_pick_best_slice <- aak1_from_het_moms_12h %>% pick_best_slices()
aak1_from_het_moms_12h_need_to_pick_best_slice$do_they_agree %>% table()

aak1_from_het_moms_24h_need_to_pick_best_slice <- aak1_from_het_moms_24h %>% pick_best_slices()
aak1_from_het_moms_24h_need_to_pick_best_slice$do_they_agree %>% table()

WT_from_het_moms_12h_need_to_pick_best_slice <- WT_from_het_moms_12h %>% pick_best_slices()
WT_from_het_moms_12h_need_to_pick_best_slice$do_they_agree %>% table()

WT_from_het_moms_24h_need_to_pick_best_slice <- WT_from_het_moms_24h %>% pick_best_slices()
WT_from_het_moms_24h_need_to_pick_best_slice$do_they_agree %>% table()


# define a function to look at cv of intensity between best slices --------

diff_best_slices <- function(df, df_need_to_pick) {
  df_need_to_pick %<>% filter(do_they_agree == FALSE)
  df %<>% filter(PGC_ID %in% df_need_to_pick$PGC_ID)
  rows_to_keep <- c()
  for(i in 1:nrow(df)) {
    if(df$slice[i] %in% as.character(filter(df_need_to_pick, PGC_ID == df$PGC_ID[i])[, c(2:3)])) {
      rows_to_keep %<>% c(i)
    }
  }
  df <- df[rows_to_keep, ]
  
  df %<>%
    group_by(PGC_ID) %>%
    summarise(diff_over_mean = diff(range(RawIntDen))/mean(RawIntDen)) %>%
    arrange(diff_over_mean)
  
  return(df)
}


# how many difference do IntDen- and area-decided slices have? ------------

aak1_12h_diff_best <- diff_best_slices(aak1_12h, aak1_12h_need_to_pick_best_slice)
aak1_24h_diff_best <- diff_best_slices(aak1_24h, aak1_24h_need_to_pick_best_slice)

WT_12h_diff_best <- diff_best_slices(WT_12h, WT_12h_need_to_pick_best_slice)
WT_24h_diff_best <- diff_best_slices(WT_24h, WT_24h_need_to_pick_best_slice)

aak1_from_het_moms_12h_diff_best <- diff_best_slices(aak1_from_het_moms_12h, aak1_from_het_moms_12h_need_to_pick_best_slice)
aak1_from_het_moms_24h_diff_best <- diff_best_slices(aak1_from_het_moms_24h, aak1_from_het_moms_24h_need_to_pick_best_slice)

WT_from_het_moms_12h_diff_best <- diff_best_slices(WT_from_het_moms_12h, WT_from_het_moms_12h_need_to_pick_best_slice)
WT_from_het_moms_24h_diff_best <- diff_best_slices(WT_from_het_moms_24h, WT_from_het_moms_24h_need_to_pick_best_slice)

aak1_12h_diff_best$diff_over_mean %>% mean()
aak1_24h_diff_best$diff_over_mean %>% mean()

WT_12h_diff_best$diff_over_mean %>% mean()
WT_24h_diff_best$diff_over_mean %>% mean()

aak1_from_het_moms_12h_diff_best$diff_over_mean %>% mean()
aak1_from_het_moms_24h_diff_best$diff_over_mean %>% mean()

WT_from_het_moms_12h_diff_best$diff_over_mean %>% mean()
WT_from_het_moms_24h_diff_best$diff_over_mean %>% mean()


# decide on the best slices --------------------------------

aak1_12h_best_slices_int <- best_slice_df(aak1_12h, aak1_12h_need_to_pick_best_slice, "intensity")
aak1_12h_best_slices_area <- best_slice_df(aak1_12h, aak1_12h_need_to_pick_best_slice, "area")

aak1_24h_best_slices_int <- best_slice_df(aak1_24h, aak1_24h_need_to_pick_best_slice, "intensity")
aak1_24h_best_slices_area <- best_slice_df(aak1_24h, aak1_24h_need_to_pick_best_slice, "area")

WT_12h_best_slices_int <- best_slice_df(WT_12h, WT_12h_need_to_pick_best_slice, "intensity")
WT_12h_best_slices_area <- best_slice_df(WT_12h, WT_12h_need_to_pick_best_slice, "area")

WT_24h_best_slices_int <- best_slice_df(WT_24h, WT_24h_need_to_pick_best_slice, "intensity")
WT_24h_best_slices_area <- best_slice_df(WT_24h, WT_24h_need_to_pick_best_slice, "area")

aak1_from_het_moms_12h_best_slices_int <- best_slice_df(aak1_from_het_moms_12h, aak1_from_het_moms_12h_need_to_pick_best_slice, "intensity")
aak1_from_het_moms_12h_best_slices_area <- best_slice_df(aak1_from_het_moms_12h, aak1_from_het_moms_12h_need_to_pick_best_slice, "area")

aak1_from_het_moms_24h_best_slices_int <- best_slice_df(aak1_from_het_moms_24h, aak1_from_het_moms_24h_need_to_pick_best_slice, "intensity")
aak1_from_het_moms_24h_best_slices_area <- best_slice_df(aak1_from_het_moms_24h, aak1_from_het_moms_24h_need_to_pick_best_slice, "area")

WT_from_het_moms_12h_best_slices_int <- best_slice_df(WT_from_het_moms_12h, WT_from_het_moms_12h_need_to_pick_best_slice, "intensity")
WT_from_het_moms_12h_best_slices_area <- best_slice_df(WT_from_het_moms_12h, WT_from_het_moms_12h_need_to_pick_best_slice, "area")

WT_from_het_moms_24h_best_slices_int <- best_slice_df(WT_from_het_moms_24h, WT_from_het_moms_24h_need_to_pick_best_slice, "intensity")
WT_from_het_moms_24h_best_slices_area <- best_slice_df(WT_from_het_moms_24h, WT_from_het_moms_24h_need_to_pick_best_slice, "area")


# Plot pilot and real expt separately -- plot proportion in the inner compartment using IntDen --------------------------------

compare_int <- data.frame(condition = c(rep("wildtype_12h_pilot", nrow(WT_12h_best_slices_int)),
                                        rep("wildtype_12h_real", nrow(WT_from_het_moms_12h_best_slices_int)),
                                        rep("wildtype_24h_pilot", nrow(WT_24h_best_slices_int)),
                                        rep("wildtype_24h_real", nrow(WT_from_het_moms_24h_best_slices_int)),
                                        rep("aak1_12h_pilot", nrow(aak1_12h_best_slices_int)),
                                        rep("aak1_12h_real", nrow(aak1_from_het_moms_12h_best_slices_int)),
                                        rep("aak1_24h_pilot", nrow(aak1_24h_best_slices_int)),
                                        rep("aak1_24h_real", nrow(aak1_from_het_moms_24h_best_slices_int))),
                          prop_Int_inner_comp = c(WT_12h_best_slices_int$prop_Int_inner_comp, 
                                                  WT_from_het_moms_12h_best_slices_int$prop_Int_inner_comp,
                                                  WT_24h_best_slices_int$prop_Int_inner_comp,
                                                  WT_from_het_moms_24h_best_slices_int$prop_Int_inner_comp,
                                                  aak1_12h_best_slices_int$prop_Int_inner_comp, 
                                                  aak1_from_het_moms_12h_best_slices_int$prop_Int_inner_comp,
                                                  aak1_24h_best_slices_int$prop_Int_inner_comp,
                                                  aak1_from_het_moms_24h_best_slices_int$prop_Int_inner_comp))

compare_int$condition %>% unique()
compare_int$condition %<>% as.factor() %>% fct_relevel(c("wildtype_12h_pilot",
                                                         "aak1_12h_pilot",
                                                         "wildtype_24h_pilot",
                                                         "aak1_24h_pilot",
                                                         "wildtype_12h_real",
                                                         "aak1_12h_real",
                                                         "wildtype_24h_real",
                                                         "aak1_24h_real"))

ggplot(data = compare_int,
       mapping = aes(x = condition,
                     y = prop_Int_inner_comp,
                     color = condition)) +
  geom_boxplot() +
  theme_classic() +
  theme(aspect.ratio = 1) +
  labs(title = "Proportion of total integrated intensity in the inner compartment ~ starvation duration\n(p-value =  using best slices defined by largest total integrated intensity)\nInner compartment is a circle whose area is 11.1% of the total ROI area\nand whose center is at the chromatin centroid")


# Plot pilot and real expt separately -- plot proportion in the inner compartment using area --------------------------------

compare_area <- data.frame(condition = c(rep("wildtype_12h_pilot", nrow(WT_12h_best_slices_area)),
                                         rep("wildtype_12h_real", nrow(WT_from_het_moms_12h_best_slices_area)),
                                         rep("wildtype_24h_pilot", nrow(WT_24h_best_slices_area)),
                                         rep("wildtype_24h_real", nrow(WT_from_het_moms_24h_best_slices_area)),
                                         rep("aak1_12h_pilot", nrow(aak1_12h_best_slices_area)),
                                         rep("aak1_12h_real", nrow(aak1_from_het_moms_12h_best_slices_area)),
                                         rep("aak1_24h_pilot", nrow(aak1_24h_best_slices_area)),
                                         rep("aak1_24h_real", nrow(aak1_from_het_moms_24h_best_slices_area))),
                           prop_Int_inner_comp = c(WT_12h_best_slices_area$prop_Int_inner_comp, 
                                                   WT_from_het_moms_12h_best_slices_area$prop_Int_inner_comp,
                                                   WT_24h_best_slices_area$prop_Int_inner_comp,
                                                   WT_from_het_moms_24h_best_slices_area$prop_Int_inner_comp,
                                                   aak1_12h_best_slices_area$prop_Int_inner_comp, 
                                                   aak1_from_het_moms_12h_best_slices_area$prop_Int_inner_comp,
                                                   aak1_24h_best_slices_area$prop_Int_inner_comp,
                                                   aak1_from_het_moms_24h_best_slices_area$prop_Int_inner_comp))

compare_area$condition %>% unique()
compare_area$condition %<>% as.factor() %>% fct_relevel(c("wildtype_12h_pilot",
                                                          "aak1_12h_pilot",
                                                          "wildtype_24h_pilot",
                                                          "aak1_24h_pilot",
                                                          "wildtype_12h_real",
                                                          "aak1_12h_real",
                                                          "wildtype_24h_real",
                                                          "aak1_24h_real"))

ggplot(data = compare_area,
       mapping = aes(x = condition,
                     y = prop_Int_inner_comp,
                     color = condition)) +
  geom_boxplot() +
  theme_classic() +
  theme(aspect.ratio = 1) +
  labs(title = "Proportion of total integrated intensity in the inner compartment ~ starvation duration\n(p-value =  using best slices defined by biggest area)\nInner compartment is a circle whose area is 11.1% of the total ROI area\nand whose center is at the chromatin centroid")


# Stats for pilot and real expt separately --------------------------------

# Bartlett's test to examine equal variances
bartlett.test(prop_Int_inner_comp ~ condition, data = compare_int) #  p-value = 0.422, can pool variance.
bartlett.test(prop_Int_inner_comp ~ condition, data = compare_area) #  p-value = 0.5921, can pool variance.

# t-test for 12h
t.test(filter(compare_int, condition == "aak1_12h_pilot")$prop_Int_inner_comp,
       filter(compare_int, condition == "wildtype_12h_pilot")$prop_Int_inner_comp,
       var.equal = TRUE) # p-value = 0.0395

t.test(filter(compare_int, condition == "aak1_12h_real")$prop_Int_inner_comp,
       filter(compare_int, condition == "wildtype_12h_real")$prop_Int_inner_comp,
       var.equal = TRUE) # p-value = 0.158

t.test(filter(compare_area, condition == "aak1_12h_pilot")$prop_Int_inner_comp,
       filter(compare_area, condition == "wildtype_12h_pilot")$prop_Int_inner_comp,
       var.equal = TRUE) # p-value = 0.04533

t.test(filter(compare_area, condition == "aak1_12h_real")$prop_Int_inner_comp,
       filter(compare_area, condition == "wildtype_12h_real")$prop_Int_inner_comp,
       var.equal = TRUE) # p-value = 0.08514

# t-test for 24h vs 12h using int
t.test(filter(compare_int, condition == "wildtype_24h_pilot")$prop_Int_inner_comp,
       filter(compare_int, condition == "wildtype_12h_pilot")$prop_Int_inner_comp,
       var.equal = TRUE) # p-value = 9.788e-06

t.test(filter(compare_int, condition == "aak1_24h_pilot")$prop_Int_inner_comp,
       filter(compare_int, condition == "aak1_12h_pilot")$prop_Int_inner_comp,
       var.equal = TRUE) # p-value = 0.5103

t.test(filter(compare_int, condition == "wildtype_24h_real")$prop_Int_inner_comp,
       filter(compare_int, condition == "wildtype_12h_real")$prop_Int_inner_comp,
       var.equal = TRUE) # p-value = 0.01693

t.test(filter(compare_int, condition == "aak1_24h_real")$prop_Int_inner_comp,
       filter(compare_int, condition == "aak1_12h_real")$prop_Int_inner_comp,
       var.equal = TRUE) # p-value = 0.3295

# t-test for 24h vs 12h using area
t.test(filter(compare_area, condition == "wildtype_24h_pilot")$prop_Int_inner_comp,
       filter(compare_area, condition == "wildtype_12h_pilot")$prop_Int_inner_comp,
       var.equal = TRUE) # p-value = 3.833e-06

t.test(filter(compare_area, condition == "aak1_24h_pilot")$prop_Int_inner_comp,
       filter(compare_area, condition == "aak1_12h_pilot")$prop_Int_inner_comp,
       var.equal = TRUE) # p-value = 0.6624

t.test(filter(compare_area, condition == "wildtype_24h_real")$prop_Int_inner_comp,
       filter(compare_area, condition == "wildtype_12h_real")$prop_Int_inner_comp,
       var.equal = TRUE) # p-value = 0.004508

t.test(filter(compare_area, condition == "aak1_24h_real")$prop_Int_inner_comp,
       filter(compare_area, condition == "aak1_12h_real")$prop_Int_inner_comp,
       var.equal = TRUE) # p-value = 0.8297


# plot proportion in the inner compartment using IntDen --------------------------------

compare_int <- data.frame(condition = c(rep("wildtype_12h_post_bleach", nrow(WT_12h_best_slices_int)),
                                        rep("wildtype_12h_post_bleach", nrow(WT_from_het_moms_12h_best_slices_int)),
                                        rep("wildtype_24h_post_bleach", nrow(WT_24h_best_slices_int)),
                                        rep("wildtype_24h_post_bleach", nrow(WT_from_het_moms_24h_best_slices_int)),
                                        rep("aak1_12h_post_bleach", nrow(aak1_12h_best_slices_int)),
                                        rep("aak1_12h_post_bleach", nrow(aak1_from_het_moms_12h_best_slices_int)),
                                        rep("aak1_24h_post_bleach", nrow(aak1_from_het_moms_24h_best_slices_int))),
                          prop_Int_inner_comp = c(WT_12h_best_slices_int$prop_Int_inner_comp, 
                                                  WT_from_het_moms_12h_best_slices_int$prop_Int_inner_comp,
                                                  WT_24h_best_slices_int$prop_Int_inner_comp,
                                                  WT_from_het_moms_24h_best_slices_int$prop_Int_inner_comp,
                                                  aak1_12h_best_slices_int$prop_Int_inner_comp, 
                                                  aak1_from_het_moms_12h_best_slices_int$prop_Int_inner_comp,
                                                  aak1_from_het_moms_24h_best_slices_int$prop_Int_inner_comp))

compare_int$condition %>% unique()
compare_int$condition %<>% as.factor() %>% fct_relevel(c("wildtype_12h_post_bleach", "aak1_12h_post_bleach", "wildtype_24h_post_bleach", "aak1_24h_post_bleach"))

ggplot(data = compare_int,
       mapping = aes(x = condition,
                     y = prop_Int_inner_comp,
                     color = condition)) +
  geom_boxplot() +
  theme_classic() +
  theme(aspect.ratio = 1) +
  labs(title = "Proportion of total integrated intensity in the inner compartment ~ starvation duration\n(p-value =  using best slices defined by largest total integrated intensity)\nInner compartment is a circle whose area is 11.1% of the total ROI area\nand whose center is at the chromatin centroid")


# plot proportion in the inner compartment using area --------------------------------

compare_area <- data.frame(condition = c(rep("wildtype_12h_post_bleach", nrow(WT_12h_best_slices_area)),
                                         rep("wildtype_12h_post_bleach", nrow(WT_from_het_moms_12h_best_slices_area)),
                                         rep("wildtype_24h_post_bleach", nrow(WT_24h_best_slices_area)),
                                         rep("wildtype_24h_post_bleach", nrow(WT_from_het_moms_24h_best_slices_area)),
                                         rep("aak1_12h_post_bleach", nrow(aak1_12h_best_slices_area)),
                                         rep("aak1_12h_post_bleach", nrow(aak1_from_het_moms_12h_best_slices_area)),
                                         rep("aak1_24h_post_bleach", nrow(aak1_from_het_moms_24h_best_slices_area))),
                           prop_Int_inner_comp = c(WT_12h_best_slices_area$prop_Int_inner_comp, 
                                                   WT_from_het_moms_12h_best_slices_area$prop_Int_inner_comp,
                                                   WT_24h_best_slices_area$prop_Int_inner_comp,
                                                   WT_from_het_moms_24h_best_slices_area$prop_Int_inner_comp,
                                                   aak1_12h_best_slices_area$prop_Int_inner_comp, 
                                                   aak1_from_het_moms_12h_best_slices_area$prop_Int_inner_comp,
                                                   aak1_from_het_moms_24h_best_slices_area$prop_Int_inner_comp))

compare_area$condition %>% unique()
compare_area$condition %<>% as.factor() %>% fct_relevel(c("wildtype_12h_post_bleach", "aak1_12h_post_bleach", "wildtype_24h_post_bleach", "aak1_24h_post_bleach"))

ggplot(data = compare_area,
       mapping = aes(x = condition,
                     y = prop_Int_inner_comp,
                     color = condition)) +
  geom_boxplot() +
  theme_classic() +
  theme(aspect.ratio = 1) +
  labs(title = "Proportion of total integrated intensity in the inner compartment ~ starvation duration\n(p-value =  using best slices defined by biggest area)\nInner compartment is a circle whose area is 11.1% of the total ROI area\nand whose center is at the chromatin centroid")


# compare proportion in the inner compartment ------------------------------------------------------------

# they all pass normality test -- can't reject the null hypothesis that they are normally distributed

shapiro.test(aak1_12h_best_slices_int$prop_Int_inner_comp) # p-value = 0.3937
shapiro.test(aak1_12h_best_slices_area$prop_Int_inner_comp) # p-value = 0.2997

shapiro.test(aak1_24h_best_slices_int$prop_Int_inner_comp) # p-value = 0.6808
shapiro.test(aak1_24h_best_slices_area$prop_Int_inner_comp) # p-value = 0.08371

shapiro.test(WT_12h_best_slices_int$prop_Int_inner_comp) # p-value = 0.7532
shapiro.test(WT_12h_best_slices_area$prop_Int_inner_comp) # p-value = 0.2755

shapiro.test(WT_24h_best_slices_int$prop_Int_inner_comp) # p-value = 0.4665
shapiro.test(WT_24h_best_slices_area$prop_Int_inner_comp) # p-value = 0.2909

shapiro.test(aak1_from_het_moms_12h_best_slices_int$prop_Int_inner_comp) # p-value = 0.5868
shapiro.test(aak1_from_het_moms_12h_best_slices_area$prop_Int_inner_comp) # p-value = 0.4971

shapiro.test(aak1_from_het_moms_24h_best_slices_int$prop_Int_inner_comp) # p-value = 0.6284
shapiro.test(aak1_from_het_moms_24h_best_slices_area$prop_Int_inner_comp) # p-value = 0.5108

shapiro.test(WT_from_het_moms_12h_best_slices_int$prop_Int_inner_comp) # p-value = 0.2253
shapiro.test(WT_from_het_moms_12h_best_slices_area$prop_Int_inner_comp) # p-value = 0.2729

shapiro.test(WT_from_het_moms_24h_best_slices_int$prop_Int_inner_comp) # p-value = 0.7398
shapiro.test(WT_from_het_moms_24h_best_slices_area$prop_Int_inner_comp) # p-value = 0.8542

# test for homogeneity of variance across groups

starvation_x_genotype_int <- compare_int

starvation_x_genotype_int$starvation <- NA
starvation_x_genotype_int$starvation[which(str_detect(starvation_x_genotype_int$condition, "12h_post_bleach"))] <- "12h_post_bleach"
starvation_x_genotype_int$starvation[which(str_detect(starvation_x_genotype_int$condition, "24h_post_bleach"))] <- "24h_post_bleach"

starvation_x_genotype_int$genotype <- NA
starvation_x_genotype_int$genotype[which(str_detect(starvation_x_genotype_int$condition, "wildtype"))] <- "wildtype"
starvation_x_genotype_int$genotype[which(str_detect(starvation_x_genotype_int$condition, "aak1"))] <- "aak1"

leveneTest(data = starvation_x_genotype_int, prop_Int_inner_comp ~ starvation*genotype) # Pr 0.08602, equal variance, can do ANOVA

starvation_x_genotype_area <- compare_area

starvation_x_genotype_area$starvation <- NA
starvation_x_genotype_area$starvation[which(str_detect(starvation_x_genotype_area$condition, "12h_post_bleach"))] <- "12h_post_bleach"
starvation_x_genotype_area$starvation[which(str_detect(starvation_x_genotype_area$condition, "24h_post_bleach"))] <- "24h_post_bleach"

starvation_x_genotype_area$genotype <- NA
starvation_x_genotype_area$genotype[which(str_detect(starvation_x_genotype_area$condition, "wildtype"))] <- "wildtype"
starvation_x_genotype_area$genotype[which(str_detect(starvation_x_genotype_area$condition, "aak1"))] <- "aak1"

leveneTest(data = starvation_x_genotype_area, prop_Int_inner_comp ~ starvation*genotype) # Pr 0.3059, equal variance, can do ANOVA

# Bartlett's test to examine equal variances
bartlett.test(prop_Int_inner_comp ~ condition, data = compare_int) #  p-value = 0.2755, can pool variance.
bartlett.test(prop_Int_inner_comp ~ condition, data = compare_area) #  p-value = 0.3398, can pool variance.

# t-test
t.test(filter(compare_int, condition == "aak1_12h_post_bleach")$prop_Int_inner_comp,
       filter(compare_int, condition == "wildtype_12h_post_bleach")$prop_Int_inner_comp,
       var.equal = TRUE) # p-value = 0.4658

t.test(filter(compare_int, condition == "aak1_24h_post_bleach")$prop_Int_inner_comp,
       filter(compare_int, condition == "wildtype_24h_post_bleach")$prop_Int_inner_comp,
       var.equal = TRUE) # p-value = 0.0445

t.test(filter(compare_int, condition == "wildtype_24h_post_bleach")$prop_Int_inner_comp,
       filter(compare_int, condition == "wildtype_12h_post_bleach")$prop_Int_inner_comp,
       var.equal = TRUE) # p-value = 9.75e-06

t.test(filter(compare_int, condition == "aak1_24h_post_bleach")$prop_Int_inner_comp,
       filter(compare_int, condition == "aak1_12h_post_bleach")$prop_Int_inner_comp,
       var.equal = TRUE) # p-value = 0.001038

t.test(filter(compare_area, condition == "aak1_12h_post_bleach")$prop_Int_inner_comp,
       filter(compare_area, condition == "wildtype_12h_post_bleach")$prop_Int_inner_comp,
       var.equal = TRUE) # p-value = 0.6142

t.test(filter(compare_area, condition == "aak1_24h_post_bleach")$prop_Int_inner_comp,
       filter(compare_area, condition == "wildtype_24h_post_bleach")$prop_Int_inner_comp,
       var.equal = TRUE) # p-value = 0.0003135

t.test(filter(compare_area, condition == "wildtype_24h_post_bleach")$prop_Int_inner_comp,
       filter(compare_area, condition == "wildtype_12h_post_bleach")$prop_Int_inner_comp,
       var.equal = TRUE) # p-value = 1.068e-06

t.test(filter(compare_area, condition == "aak1_24h_post_bleach")$prop_Int_inner_comp,
       filter(compare_area, condition == "aak1_12h_post_bleach")$prop_Int_inner_comp,
       var.equal = TRUE) # p-value = 0.07587

# two-way ANOVA
aov(data = starvation_x_genotype_int, prop_Int_inner_comp ~ starvation*genotype) %>% summary() # starvation:genotype interaction 0.3635
aov(data = starvation_x_genotype_area, prop_Int_inner_comp ~ starvation*genotype) %>% summary() # starvation:genotype interaction 0.02196

starvation_x_genotype_area %>% dplyr::group_by(condition) %>% dplyr::summarise(num_obs = n())

# one-way ANOVA + post hoc Tukey
aov(data = compare_int, formula = prop_Int_inner_comp ~ condition) %>% TukeyHSD(., which = "condition")
# wildtype_24h_post_bleach-wildtype_12h_post_bleach p adj 0.0000121
# wildtype_24h_post_bleach-aak1_12h_post_bleach p adj 0.0000001
# aak1_24h_post_bleach-aak1_12h_post_bleach p adj 0.0132290

aov(data = compare_area, formula = prop_Int_inner_comp ~ condition) %>% TukeyHSD(., which = "condition")
# wildtype_24h_post_bleach-wildtype_12h_post_bleach p adj 0.0000010
# wildtype_24h_post_bleach-aak1_12h_post_bleach p adj 0.0000000
# aak1_24h_post_bleach-wildtype_24h_post_bleach p adj 0.0010635


# define a function to combine separate line profiling CSVs --------------------------------

combine_profiling_csv <- function(rootdir) {
  # set working directory
  setwd("C:/Users/jingx/Duke Bio_Ea Dropbox/Baugh Lab/2_baugh_lab_users/Jingxian/confocal_Zeiss880")
  setwd(rootdir)
  
  # get all intensity line profiling CSV files
  # full.names: the path (relative to the current working directory) is prepended to the file names
  # recursive: go into subdirectories
  # include.dirs: don't just return files. Also return directories.
  line_profiling_file_names <- 
    list.files(pattern = "- Copy - Copy", full.names = TRUE, recursive = FALSE, include.dirs = TRUE) %>% 
    list.files(pattern = "pgc", full.names = TRUE, recursive = FALSE, include.dirs = TRUE) %>%
    list.files(pattern = "_intensity_profile.csv", full.names = TRUE, recursive = TRUE, include.dirs = FALSE) %>%
    sapply(FUN = function(x) {unlist(str_remove_all(x, "\\.\\/"))})
  
  print(length(line_profiling_file_names)) %>% paste(., "names have been retrieved")
  
  # create a dataframe where you will save all the passed worms' values
  line_profiling <- data.frame(raw_int = NA,
                               nrow = NA,
                               dir = NA,
                               int_norm_by_max = NA,
                               prop_of_total_scan_distance = NA,
                               PGC_ID = NA,
                               slice = NA)
  
  # loop through all the CSV files
  for (j in 1:length(line_profiling_file_names)) {
    df <- read_csv(file = line_profiling_file_names[j], col_names = TRUE) %>% .[, -1] %>% unique()
    df$dir <- line_profiling_file_names[j]
    max_int <- df$raw_int %>% max()
    df$int_norm_by_max <- df$raw_int/max_int
    df$prop_of_total_scan_distance <- df$nrow/tail(df$nrow, 1)
    
    df$PGC_ID <- NA
    df$slice <- NA
    
    for (i in 1:nrow(df)) {
      temp <- df$dir[i] %>% str_split("/") %>% unlist()
      image <- temp[length(temp) - 2]
      image %<>% str_split("_Out") %>% unlist() %>% .[1]
      pgc <- temp[length(temp) - 1]
      df$PGC_ID[i] <- paste(image, pgc, sep = " ")
      
      slice_id <- temp[length(temp)]
      slice_id %<>% str_split("Out-|Image-|.tif") %>% unlist() %>% .[2]
      df$slice[i] <- slice_id
    }
    
    line_profiling %<>% rbind(df)
  }
  
  # should get row number minus 1, because the first row is all-NA
  line_profiling %<>%
    na.omit() %>%
    select(c("PGC_ID", "slice", "prop_of_total_scan_distance", "int_norm_by_max"))
  
  return(line_profiling)
}


# aggregate line intensity profiling results ------------------------------------

aak1_12h_profiling <- combine_profiling_csv(rootdir = "./20241120_aak-1_GFP_neg_12h_chr_compaction/processed_auto_3D")
aak1_24h_profiling <- combine_profiling_csv(rootdir = "./20241120_aak-1_GFP_neg_24h_chr_compaction/processed_auto_3D")

WT_12h_profiling <- combine_profiling_csv(rootdir = "./20241120_WT_GFP_neg_12h_chr_compaction/processed_auto_3D")
WT_24h_profiling <- combine_profiling_csv(rootdir = "./20241120_WT_GFP_neg_24h_chr_compaction/processed_auto_3D")

aak1_from_het_moms_12h_profiling <- combine_profiling_csv(rootdir = "./20241127_aak-1_GFP_neg_12h_chr_compaction/processed_auto_3D")
aak1_from_het_moms_24h_profiling <- combine_profiling_csv(rootdir = "./20241127_aak-1_GFP_neg_24h_chr_compaction/processed_auto_3D")

WT_from_het_moms_12h_profiling <- combine_profiling_csv(rootdir = "./20241127_WT_GFP_neg_12h_chr_compaction/processed_auto_3D")
WT_from_het_moms_24h_profiling <- combine_profiling_csv(rootdir = "./20241127_WT_GFP_neg_24h_chr_compaction/processed_auto_3D")


# define a function to extract best slices' line profiling results -----------

best_slices_line_profiling <- function(df, best_slices_df, conf = 0.95) {
  
  df %<>% mutate(id = paste(PGC_ID, slice, sep = "_"))
  best_slices_df %<>% mutate(id = paste(PGC_ID, slice, sep = "_"))
  
  df %<>% filter(id %in% best_slices_df$id)
  
  df %<>%
    mutate(prop_of_total_scan_distance = round(prop_of_total_scan_distance, 1)) %>%
    group_by(prop_of_total_scan_distance) %>%
    reframe(int_norm_by_max_mean = mean(int_norm_by_max),
            lower_CI = t.test(int_norm_by_max, conf.level = conf)$conf.int[1],
            upper_CI = t.test(int_norm_by_max, conf.level = conf)$conf.int[2])
  
  colnames(df)[2] <- "int_norm_by_max"
  
  # df$lower_CI[is.na(df$lower_CI)] <- df$int_norm_by_max[is.na(df$lower_CI)]
  # df$upper_CI[is.na(df$upper_CI)] <- df$int_norm_by_max[is.na(df$upper_CI)]
  # 
  # df$lower_CI[df$lower_CI < 0] <- 0
  # df$lower_CI[df$lower_CI > 1] <- 1
  # 
  # df$upper_CI[df$upper_CI < 0] <- 0
  # df$upper_CI[df$upper_CI > 1] <- 1
  
  return(df)

}


# extract best slices' line profiling results -----------------------------

aak1_12h_profiling_best_slices_int <- best_slices_line_profiling(aak1_12h_profiling, aak1_12h_best_slices_int, conf = 0.95)
aak1_24h_profiling_best_slices_int <- best_slices_line_profiling(aak1_24h_profiling, aak1_24h_best_slices_int, conf = 0.95)

WT_12h_profiling_best_slices_int <- best_slices_line_profiling(WT_12h_profiling, WT_12h_best_slices_int, conf = 0.95)
WT_24h_profiling_best_slices_int <- best_slices_line_profiling(WT_24h_profiling, WT_24h_best_slices_int, conf = 0.95)

aak1_from_het_moms_12h_profiling_best_slices_int <- best_slices_line_profiling(aak1_from_het_moms_12h_profiling, aak1_from_het_moms_12h_best_slices_int, conf = 0.95)
aak1_from_het_moms_24h_profiling_best_slices_int <- best_slices_line_profiling(aak1_from_het_moms_24h_profiling, aak1_from_het_moms_24h_best_slices_int, conf = 0.95)

WT_from_het_moms_12h_profiling_best_slices_int <- best_slices_line_profiling(WT_from_het_moms_12h_profiling, WT_from_het_moms_12h_best_slices_int, conf = 0.95)
WT_from_het_moms_24h_profiling_best_slices_int <- best_slices_line_profiling(WT_from_het_moms_24h_profiling, WT_from_het_moms_24h_best_slices_int, conf = 0.95)

aak1_12h_profiling_best_slices_area <- best_slices_line_profiling(aak1_12h_profiling, aak1_12h_best_slices_area, conf = 0.95)
aak1_24h_profiling_best_slices_area <- best_slices_line_profiling(aak1_24h_profiling, aak1_24h_best_slices_area, conf = 0.95)

WT_12h_profiling_best_slices_area <- best_slices_line_profiling(WT_12h_profiling, WT_12h_best_slices_area, conf = 0.95)
WT_24h_profiling_best_slices_area <- best_slices_line_profiling(WT_24h_profiling, WT_24h_best_slices_area, conf = 0.95)

aak1_from_het_moms_12h_profiling_best_slices_area <- best_slices_line_profiling(aak1_from_het_moms_12h_profiling, aak1_from_het_moms_12h_best_slices_area, conf = 0.95)
aak1_from_het_moms_24h_profiling_best_slices_area <- best_slices_line_profiling(aak1_from_het_moms_24h_profiling, aak1_from_het_moms_24h_best_slices_area, conf = 0.95)

WT_from_het_moms_12h_profiling_best_slices_area <- best_slices_line_profiling(WT_from_het_moms_12h_profiling, WT_from_het_moms_12h_best_slices_area, conf = 0.95)
WT_from_het_moms_24h_profiling_best_slices_area <- best_slices_line_profiling(WT_from_het_moms_24h_profiling, WT_from_het_moms_24h_best_slices_area, conf = 0.95)


# Plot pilot and real expt separately -- plot line profiling results ---------------------------------------------

line_profiling_compare_int <- rbind(WT_12h_profiling_best_slices_int,
                                    WT_from_het_moms_12h_profiling_best_slices_int,
                                    WT_24h_profiling_best_slices_int,
                                    WT_from_het_moms_24h_profiling_best_slices_int,
                                    aak1_12h_profiling_best_slices_int,
                                    aak1_from_het_moms_12h_profiling_best_slices_int,
                                    aak1_24h_profiling_best_slices_int,
                                    aak1_from_het_moms_24h_profiling_best_slices_int)

line_profiling_compare_int$condition <- c(rep("wildtype_12h_pilot", nrow(WT_12h_profiling_best_slices_int)),
                                          rep("wildtype_12h_real", nrow(WT_from_het_moms_12h_profiling_best_slices_int)),
                                          rep("wildtype_24h_pilot", nrow(WT_24h_profiling_best_slices_int)),
                                          rep("wildtype_24h_real", nrow(WT_from_het_moms_24h_profiling_best_slices_int)),
                                          rep("aak1_12h_pilot", nrow(aak1_12h_profiling_best_slices_int)),
                                          rep("aak1_12h_real", nrow(aak1_from_het_moms_12h_profiling_best_slices_int)),
                                          rep("aak1_24h_pilot", nrow(aak1_24h_profiling_best_slices_int)),
                                          rep("aak1_24h_real", nrow(aak1_from_het_moms_24h_profiling_best_slices_int)))

line_profiling_compare_int$condition %>% unique()
line_profiling_compare_int$condition %<>% as.factor() %>% fct_relevel(c("wildtype_12h_pilot",
                                                                        "aak1_12h_pilot",
                                                                        "wildtype_24h_pilot",
                                                                        "aak1_24h_pilot",
                                                                        "wildtype_12h_real",
                                                                        "aak1_12h_real",
                                                                        "wildtype_24h_real",
                                                                        "aak1_24h_real"))

line_profiling_compare_area <- rbind(WT_12h_profiling_best_slices_area,
                                     WT_from_het_moms_12h_profiling_best_slices_area,
                                     WT_24h_profiling_best_slices_area,
                                     WT_from_het_moms_24h_profiling_best_slices_area,
                                     aak1_12h_profiling_best_slices_area,
                                     aak1_from_het_moms_12h_profiling_best_slices_area,
                                     aak1_24h_profiling_best_slices_area,
                                     aak1_from_het_moms_24h_profiling_best_slices_area)

line_profiling_compare_area$condition <- c(rep("wildtype_12h_pilot", nrow(WT_12h_profiling_best_slices_area)),
                                           rep("wildtype_12h_real", nrow(WT_from_het_moms_12h_profiling_best_slices_area)),
                                           rep("wildtype_24h_pilot", nrow(WT_24h_profiling_best_slices_area)),
                                           rep("wildtype_24h_real", nrow(WT_from_het_moms_24h_profiling_best_slices_area)),
                                           rep("aak1_12h_pilot", nrow(aak1_12h_profiling_best_slices_area)),
                                           rep("aak1_12h_real", nrow(aak1_from_het_moms_12h_profiling_best_slices_area)),
                                           rep("aak1_24h_pilot", nrow(aak1_24h_profiling_best_slices_area)),
                                           rep("aak1_24h_real", nrow(aak1_from_het_moms_24h_profiling_best_slices_area)))

line_profiling_compare_area$condition %>% unique()
line_profiling_compare_area$condition %<>% as.factor() %>% fct_relevel(c("wildtype_12h_pilot",
                                                                         "aak1_12h_pilot",
                                                                         "wildtype_24h_pilot",
                                                                         "aak1_24h_pilot",
                                                                         "wildtype_12h_real",
                                                                         "aak1_12h_real",
                                                                         "wildtype_24h_real",
                                                                         "aak1_24h_real"))

# plot
conf = 0.95

ggplot(data = filter(line_profiling_compare_area, condition %in% c("wildtype_24h_pilot", "wildtype_12h_pilot")),
       mapping = aes(x = prop_of_total_scan_distance,
                     y = int_norm_by_max)) +
  geom_ribbon(mapping = aes(ymin = lower_CI, ymax = upper_CI, fill = condition), alpha = 0.25, color = NA) +
  geom_line(mapping = aes(color = condition), linewidth = 1) +
  theme_classic() +
  scale_x_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, 0.1),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, 0.1),
                     expand = expansion(mult = c(0, 0))) +
  theme(aspect.ratio = 1) +
  labs(title = paste("line is average trend, shade is CI", conf, sep = " "))

ggplot(data = filter(line_profiling_compare_area, condition %in% c("wildtype_24h_real", "wildtype_12h_real")),
       mapping = aes(x = prop_of_total_scan_distance,
                     y = int_norm_by_max)) +
  geom_ribbon(mapping = aes(ymin = lower_CI, ymax = upper_CI, fill = condition), alpha = 0.25, color = NA) +
  geom_line(mapping = aes(color = condition), linewidth = 1) +
  theme_classic() +
  scale_x_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, 0.1),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, 0.1),
                     expand = expansion(mult = c(0, 0))) +
  theme(aspect.ratio = 1) +
  labs(title = paste("line is average trend, shade is CI", conf, sep = " "))

ggplot(data = filter(line_profiling_compare_area, condition %in% c("aak1_24h_pilot", "aak1_12h_pilot")),
       mapping = aes(x = prop_of_total_scan_distance,
                     y = int_norm_by_max)) +
  geom_ribbon(mapping = aes(ymin = lower_CI, ymax = upper_CI, fill = condition), alpha = 0.25, color = NA) +
  geom_line(mapping = aes(color = condition), linewidth = 1) +
  theme_classic() +
  scale_x_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, 0.1),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, 0.1),
                     expand = expansion(mult = c(0, 0))) +
  theme(aspect.ratio = 1) +
  labs(title = paste("line is average trend, shade is CI", conf, sep = " "))

ggplot(data = filter(line_profiling_compare_area, condition %in% c("aak1_24h_real", "aak1_12h_real")),
       mapping = aes(x = prop_of_total_scan_distance,
                     y = int_norm_by_max)) +
  geom_ribbon(mapping = aes(ymin = lower_CI, ymax = upper_CI, fill = condition), alpha = 0.25, color = NA) +
  geom_line(mapping = aes(color = condition), linewidth = 1) +
  theme_classic() +
  scale_x_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, 0.1),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, 0.1),
                     expand = expansion(mult = c(0, 0))) +
  theme(aspect.ratio = 1) +
  labs(title = paste("line is average trend, shade is CI", conf, sep = " "))


# plot line profiling results ---------------------------------------------

line_profiling_compare_int <- rbind(WT_12h_profiling_best_slices_int,
                                    WT_from_het_moms_12h_profiling_best_slices_int,
                                    WT_24h_profiling_best_slices_int,
                                    WT_from_het_moms_24h_profiling_best_slices_int,
                                    aak1_12h_profiling_best_slices_int,
                                    aak1_from_het_moms_12h_profiling_best_slices_int,
                                    aak1_from_het_moms_24h_profiling_best_slices_int)

line_profiling_compare_int$condition <- c(rep("wildtype_12h_post_bleach", nrow(WT_12h_profiling_best_slices_int)),
                                          rep("wildtype_12h_post_bleach", nrow(WT_from_het_moms_12h_profiling_best_slices_int)),
                                          rep("wildtype_24h_post_bleach", nrow(WT_24h_profiling_best_slices_int)),
                                          rep("wildtype_24h_post_bleach", nrow(WT_from_het_moms_24h_profiling_best_slices_int)),
                                          rep("aak1_12h_post_bleach", nrow(aak1_12h_profiling_best_slices_int)),
                                          rep("aak1_12h_post_bleach", nrow(aak1_from_het_moms_12h_profiling_best_slices_int)),
                                          rep("aak1_24h_post_bleach", nrow(aak1_from_het_moms_24h_profiling_best_slices_int)))

# line_profiling_compare_int <- rbind(WT_from_het_moms_12h_profiling_best_slices_int,
#                                     WT_from_het_moms_24h_profiling_best_slices_int,
#                                     aak1_from_het_moms_12h_profiling_best_slices_int,
#                                     aak1_from_het_moms_24h_profiling_best_slices_int)
# 
# line_profiling_compare_int$condition <- c(rep("wildtype_12h_post_bleach", nrow(WT_from_het_moms_12h_profiling_best_slices_int)),
#                                           rep("wildtype_24h_post_bleach", nrow(WT_from_het_moms_24h_profiling_best_slices_int)),
#                                           rep("aak1_12h_post_bleach", nrow(aak1_from_het_moms_12h_profiling_best_slices_int)),
#                                           rep("aak1_24h_post_bleach", nrow(aak1_from_het_moms_24h_profiling_best_slices_int)))

line_profiling_compare_int$condition %>% unique()
line_profiling_compare_int$condition %<>% as.factor() %>% fct_relevel(c("wildtype_12h_post_bleach", "aak1_12h_post_bleach", "wildtype_24h_post_bleach", "aak1_24h_post_bleach"))

line_profiling_compare_area <- rbind(WT_12h_profiling_best_slices_area,
                                     WT_from_het_moms_12h_profiling_best_slices_area,
                                     WT_24h_profiling_best_slices_area,
                                     WT_from_het_moms_24h_profiling_best_slices_area,
                                     aak1_12h_profiling_best_slices_area,
                                     aak1_from_het_moms_12h_profiling_best_slices_area,
                                     aak1_from_het_moms_24h_profiling_best_slices_area)

line_profiling_compare_area$condition <- c(rep("wildtype_12h_post_bleach", nrow(WT_12h_profiling_best_slices_area)),
                                           rep("wildtype_12h_post_bleach", nrow(WT_from_het_moms_12h_profiling_best_slices_area)),
                                           rep("wildtype_24h_post_bleach", nrow(WT_24h_profiling_best_slices_area)),
                                           rep("wildtype_24h_post_bleach", nrow(WT_from_het_moms_24h_profiling_best_slices_area)),
                                           rep("aak1_12h_post_bleach", nrow(aak1_12h_profiling_best_slices_area)),
                                           rep("aak1_12h_post_bleach", nrow(aak1_from_het_moms_12h_profiling_best_slices_area)),
                                           rep("aak1_24h_post_bleach", nrow(aak1_from_het_moms_24h_profiling_best_slices_area)))

# line_profiling_compare_area <- rbind(WT_from_het_moms_12h_profiling_best_slices_area,
#                                      WT_from_het_moms_24h_profiling_best_slices_area,
#                                      aak1_from_het_moms_12h_profiling_best_slices_area,
#                                      aak1_from_het_moms_24h_profiling_best_slices_area)
# 
# line_profiling_compare_area$condition <- c(rep("wildtype_12h_post_bleach", nrow(WT_from_het_moms_12h_profiling_best_slices_area)),
#                                            rep("wildtype_24h_post_bleach", nrow(WT_from_het_moms_24h_profiling_best_slices_area)),
#                                            rep("aak1_12h_post_bleach", nrow(aak1_from_het_moms_12h_profiling_best_slices_area)),
#                                            rep("aak1_24h_post_bleach", nrow(aak1_from_het_moms_24h_profiling_best_slices_area)))
line_profiling_compare_area$condition %>% unique()
line_profiling_compare_area$condition %<>% as.factor() %>% fct_relevel(c("wildtype_12h_post_bleach", "aak1_12h_post_bleach", "wildtype_24h_post_bleach", "aak1_24h_post_bleach"))


cairo_pdf(filename = "line_profiling.pdf", width = 10, height = 10, onefile = TRUE)
# plot WT 24h post bleach vs 12 h
p1 <-
  ggplot(data = filter(line_profiling_compare_area, condition %in% c("wildtype_24h_post_bleach", "wildtype_12h_post_bleach")),
       mapping = aes(x = prop_of_total_scan_distance,
                     y = int_norm_by_max)) +
  geom_ribbon(mapping = aes(ymin = lower_CI, ymax = upper_CI, fill = condition), alpha = 0.25, color = NA) +
  geom_line(mapping = aes(color = condition), linewidth = 1) +
  theme_classic() +
  scale_x_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, 0.1),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, 0.1),
                     expand = expansion(mult = c(0, 0))) +
  theme(aspect.ratio = 1) +
  labs(title = paste("line is average trend, shade is CI", conf, sep = " "))

print(p1)

# plot aak-1 24h post bleach vs 12 h
p2 <-
  ggplot(data = filter(line_profiling_compare_area, condition %in% c("aak1_24h_post_bleach", "aak1_12h_post_bleach")),
       mapping = aes(x = prop_of_total_scan_distance,
                     y = int_norm_by_max)) +
  geom_ribbon(mapping = aes(ymin = lower_CI, ymax = upper_CI, fill = condition), alpha = 0.25, color = NA) +
  geom_line(mapping = aes(color = condition), linewidth = 1) +
  theme_classic() +
  scale_x_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, 0.1),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, 0.1),
                     expand = expansion(mult = c(0, 0))) +
  theme(aspect.ratio = 1) +
  labs(title = paste("line is average trend, shade is CI", conf, sep = " "))

print(p2)

# plot aak-1 vs WT at 12h post bleach
p3<-
  ggplot(data = filter(line_profiling_compare_area, condition %in% c("aak1_12h_post_bleach", "wildtype_12h_post_bleach")),
       mapping = aes(x = prop_of_total_scan_distance,
                     y = int_norm_by_max)) +
  geom_ribbon(mapping = aes(ymin = lower_CI, ymax = upper_CI, fill = condition), alpha = 0.25, color = NA) +
  geom_line(mapping = aes(color = condition), linewidth = 1) +
  theme_classic() +
  scale_x_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, 0.1),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, 0.1),
                     expand = expansion(mult = c(0, 0))) +
  theme(aspect.ratio = 1) +
  labs(title = paste("line is average trend, shade is CI", conf, sep = " "))

print(p3)

# plot aak-1 vs WT at 24h post bleach
p4 <-
  ggplot(data = filter(line_profiling_compare_area, condition %in% c("aak1_24h_post_bleach", "wildtype_24h_post_bleach")),
       mapping = aes(x = prop_of_total_scan_distance,
                     y = int_norm_by_max)) +
  geom_ribbon(mapping = aes(ymin = lower_CI, ymax = upper_CI, fill = condition), alpha = 0.25, color = NA) +
  geom_line(mapping = aes(color = condition), linewidth = 1) +
  theme_classic() +
  scale_x_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, 0.1),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, 0.1),
                     expand = expansion(mult = c(0, 0))) +
  theme(aspect.ratio = 1) +
  labs(title = paste("line is average trend, shade is CI", conf, sep = " "))

print(p4)

dev.off()


# for pilot, plot proportion in the inner compartment using IntDen --------------------------------

pilot_compare_int <- data.frame(condition = c(rep("wildtype_12h_post_bleach", nrow(WT_12h_best_slices_int)),
                                              rep("wildtype_24h_post_bleach", nrow(WT_24h_best_slices_int)),
                                              rep("aak1_12h_post_bleach", nrow(aak1_12h_best_slices_int)),
                                              rep("aak1_24h_post_bleach", nrow(aak1_24h_best_slices_int))),
                                prop_Int_inner_comp = c(WT_12h_best_slices_int$prop_Int_inner_comp, 
                                                        WT_24h_best_slices_int$prop_Int_inner_comp,
                                                        aak1_12h_best_slices_int$prop_Int_inner_comp, 
                                                        aak1_24h_best_slices_int$prop_Int_inner_comp))

pilot_compare_int$condition %>% unique()
pilot_compare_int$condition %<>% as.factor() %>% fct_relevel(c("wildtype_12h_post_bleach", "aak1_12h_post_bleach", "wildtype_24h_post_bleach", "aak1_24h_post_bleach"))

ggplot(data = pilot_compare_int,
       mapping = aes(x = condition,
                     y = prop_Int_inner_comp,
                     color = condition)) +
  geom_boxplot() +
  theme_classic() +
  theme(aspect.ratio = 1) +
  labs(title = "Proportion of total integrated intensity in the inner compartment ~ starvation duration\n(p-value =  using best slices defined by largest total integrated intensity)\nInner compartment is a circle whose area is 11.1% of the total ROI area\nand whose center is at the chromatin centroid")


# for pilot, plot proportion in the inner compartment using area --------------------------------

pilot_compare_area <- data.frame(condition = c(rep("wildtype_12h_post_bleach", nrow(WT_12h_best_slices_area)),
                                               rep("wildtype_24h_post_bleach", nrow(WT_24h_best_slices_area)),
                                               rep("aak1_12h_post_bleach", nrow(aak1_12h_best_slices_area)),
                                               rep("aak1_24h_post_bleach", nrow(aak1_24h_best_slices_area))),
                                 prop_Int_inner_comp = c(WT_12h_best_slices_area$prop_Int_inner_comp, 
                                                         WT_24h_best_slices_area$prop_Int_inner_comp,
                                                         aak1_12h_best_slices_area$prop_Int_inner_comp, 
                                                         aak1_24h_best_slices_area$prop_Int_inner_comp))

pilot_compare_area$condition %>% unique()
pilot_compare_area$condition %<>% as.factor() %>% fct_relevel(c("wildtype_12h_post_bleach", "aak1_12h_post_bleach", "wildtype_24h_post_bleach", "aak1_24h_post_bleach"))

ggplot(data = pilot_compare_area,
       mapping = aes(x = condition,
                     y = prop_Int_inner_comp,
                     color = condition)) +
  geom_boxplot() +
  theme_classic() +
  theme(aspect.ratio = 1) +
  labs(title = "Proportion of total integrated intensity in the inner compartment ~ starvation duration\n(p-value =  using best slices defined by biggest area)\nInner compartment is a circle whose area is 11.1% of the total ROI area\nand whose center is at the chromatin centroid")





