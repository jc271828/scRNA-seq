library(tidyverse)
library(magrittr)
library(cowplot)

setwd("../determ_mt_prop_cutoff/")

ws286_BestExt_files <- list.files()[list.files() %>% str_detect("\\.csv")]

for (i in 1:length(ws286_BestExt_files)){
  assign(str_remove_all(ws286_BestExt_files[i], "\\.csv"), read_csv(ws286_BestExt_files[i], col_names = T))
}

# plot CDF of UMI count distribution of all barcodes
setwd("../determ_minUMI/")

# combine samples and plot
ws286_BestExt <- rbind(fed0_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc,
                       fed1_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc,
                       fed2_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc,
                       fed3_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc,
                       `fed4-1_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc`,
                       `fed4-2_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc`,
                       fed5_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc,
                       `fed6-1_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc`,
                       `fed6-2_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc`,
                       starv0_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc,
                       starv1_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc,
                       starv2_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc,
                       starv3_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc,
                       `starv4-2_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc`,
                       starv5_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc,
                       `starv6-1_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc`,
                       `starv6-2_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc`)

# pick min.UMI.low.quality
cairo_pdf(onefile = T, filename = "UmiDistrWs286BestExt20230422.pdf", width = 10, height = 10)

p_ws286_BestExt <- 
  ggplot(ws286_BestExt, aes(total)) +
  stat_ecdf() +
  theme_classic() +
  theme(aspect.ratio = 1) +
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0)),
                     breaks = seq(0, 1, 0.05)) +
  labs(title = "UMI count per barcode of all samples combined, mapped to best ext",
       x = "UMI count of a barcode")

p_ws286_BestExt
p_ws286_BestExt + coord_cartesian(xlim = c(0, 10000))
p_ws286_BestExt + coord_cartesian(xlim = c(0, 2500)) + scale_x_continuous(breaks = seq(0, 2500, 100), expand = expansion(mult = c(0, 0)))
p_ws286_BestExt + coord_cartesian(xlim = c(0, 1000)) + scale_x_continuous(breaks = seq(0, 1000, 50), expand = expansion(mult = c(0, 0)))
p_ws286_BestExt + coord_cartesian(xlim = c(0, 250)) + scale_x_continuous(breaks = seq(0, 250, 10), expand = expansion(mult = c(0, 0)))
p_ws286_BestExt + coord_cartesian(xlim = c(0, 50)) + scale_x_continuous(breaks = seq(0, 50, 1), expand = expansion(mult = c(0, 0)))

plot_count_distr <- function(df) {
  str <- deparse(substitute(df)) %>% str_split("_") %>% unlist() %>% .[1]
  p <-
    ggplot(df, aes(total)) +
    stat_ecdf() +
    theme_classic() +
    theme(aspect.ratio = 1) +
    coord_cartesian(xlim = c(0, 50)) +
    scale_x_continuous(expand = expansion(mult = c(0, 0))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0)),
                       breaks = seq(from = 0, to = 1, by = 0.05)) +
    labs(title = str,
         x = "UMI count of a barcode")
  
  return(p)
}

plot_grid(plot_count_distr(fed0_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc),
          plot_count_distr(fed1_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc),
          plot_count_distr(fed2_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc),
          plot_count_distr(fed3_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc),
          plot_count_distr(`fed4-1_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc`),
          plot_count_distr(`fed4-2_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc`),
          plot_count_distr(fed5_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc),
          plot_count_distr(`fed6-1_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc`),
          plot_count_distr(`fed6-2_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc`),
          # labels = "AUTO",
          # align = "hv",
          nrow = 3,
          ncol = 3)

plot_grid(plot_count_distr(starv0_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc),
          plot_count_distr(starv1_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc),
          plot_count_distr(starv2_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc),
          plot_count_distr(starv3_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc),
          plot_count_distr(`starv4-2_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc`),
          plot_count_distr(starv5_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc),
          plot_count_distr(`starv6-1_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc`),
          plot_count_distr(`starv6-2_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc`),
          # labels = "AUTO",
          # align = "hv",
          nrow = 3,
          ncol = 3)

dev.off()
# pick min.UMI.low.quality at 10

# pick min.UMI.high.quality
cairo_pdf(onefile = T, height = 10, width = 10, filename = "UmiAtLeast10DistrWs286BestExt20230422.pdf")

p_ws286_BestExt_UmiAtLeast10 <- 
  ggplot(filter(ws286_BestExt, total >= 10), aes(total)) +
  stat_ecdf() +
  theme_classic() +
  theme(aspect.ratio = 1) +
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0)),
                     breaks = seq(0, 1, 0.05)) +
  labs(title = "UMI count per barcode of all samples combined, mapped to best ext\nrestrictions: UMI >= 10",
       x = "UMI count of a barcode")

p_ws286_BestExt_UmiAtLeast10 
p_ws286_BestExt_UmiAtLeast10 + coord_cartesian(xlim = c(0, 10000))
p_ws286_BestExt_UmiAtLeast10 + coord_cartesian(xlim = c(0, 2500)) + scale_x_continuous(breaks = seq(0, 2500, 100), expand = expansion(mult = c(0, 0)))
p_ws286_BestExt_UmiAtLeast10 + coord_cartesian(xlim = c(0, 1000)) + scale_x_continuous(breaks = seq(0, 1000, 50), expand = expansion(mult = c(0, 0)))
p_ws286_BestExt_UmiAtLeast10 + coord_cartesian(xlim = c(0, 250)) + scale_x_continuous(breaks = seq(0, 250, 10), expand = expansion(mult = c(0, 0)))
p_ws286_BestExt_UmiAtLeast10 + coord_cartesian(xlim = c(0, 50), ylim = c(0, 0.3)) + scale_x_continuous(breaks = seq(0, 50, 1), expand = expansion(mult = c(0, 0)))

plot_count_distr_AtLeast10 <- function(df) {
  str <- deparse(substitute(df)) %>% str_split("_") %>% unlist() %>% .[1]
  p <-
    ggplot(filter(df, total > 10), aes(total)) +
    stat_ecdf() +
    theme_classic() +
    theme(aspect.ratio = 1,
          axis.text.x = element_text(size = 7, angle = 20, vjust = 1, hjust=1)) +
    coord_cartesian(xlim = c(0, 1000)) +
    scale_x_continuous(expand = expansion(mult = c(0, 0)),
                       breaks = seq(0, 1000, 50)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0)),
                       breaks = seq(from = 0, to = 1, by = 0.05)) +
    labs(title = str,
         x = "UMI count of a barcode")
  
  return(p)
}

plot_grid(plot_count_distr_AtLeast10(fed0_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc),
          plot_count_distr_AtLeast10(fed1_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc),
          plot_count_distr_AtLeast10(fed2_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc),
          plot_count_distr_AtLeast10(fed3_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc),
          plot_count_distr_AtLeast10(`fed4-1_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc`),
          plot_count_distr_AtLeast10(`fed4-2_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc`),
          plot_count_distr_AtLeast10(fed5_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc),
          plot_count_distr_AtLeast10(`fed6-1_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc`),
          plot_count_distr_AtLeast10(`fed6-2_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc`),
          # labels = "AUTO",
          # align = "hv",
          nrow = 3,
          ncol = 3)

plot_grid(plot_count_distr_AtLeast10(starv0_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc),
          plot_count_distr_AtLeast10(starv1_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc),
          plot_count_distr_AtLeast10(starv2_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc),
          plot_count_distr_AtLeast10(starv3_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc),
          plot_count_distr_AtLeast10(`starv4-2_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc`),
          plot_count_distr_AtLeast10(starv5_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc),
          plot_count_distr_AtLeast10(`starv6-1_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc`),
          plot_count_distr_AtLeast10(`starv6-2_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc`),
          # labels = "AUTO",
          # align = "hv",
          nrow = 3,
          ncol = 3)

dev.off()
# pick 25 as min.UMI.high.quality in COPILOT



