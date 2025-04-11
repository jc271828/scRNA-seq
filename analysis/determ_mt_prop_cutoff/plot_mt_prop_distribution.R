library(tidyverse)
library(magrittr)
library(cowplot)

ws286_BestExt_files <- list.files()[list.files() %>% str_detect("\\.csv")]

for (i in 1:length(ws286_BestExt_files)){
  assign(str_remove_all(ws286_BestExt_files[i], "\\.csv"), read_csv(ws286_BestExt_files[i], col_names = T))
}

# plot CDF of mt_genes percentage
cairo_pdf(onefile = T, filename = "MtPropWs286BestExt20230422.pdf", width = 10, height = 10)

plot_mt_genes <- function(df) {
  str <- deparse(substitute(df)) %>% str_split("_") %>% unlist() %>% .[1]
  p <-
    ggplot(df, aes(mt_prop)) +
    stat_ecdf() +
    theme_classic() +
    theme(aspect.ratio = 1) +
    # coord_cartesian(xlim = c(0, 0.5)) +
    scale_x_continuous(expand = expansion(mult = c(0, 0)),
                       breaks = seq(from = 0, to = 1, by = 0.1)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0)),
                       breaks = seq(from = 0, to = 1, by = 0.05)) +
    labs(title = str)
    
  return(p)
}

plot_grid(plot_mt_genes(fed0_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc),
          plot_mt_genes(fed1_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc),
          plot_mt_genes(fed2_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc),
          plot_mt_genes(fed3_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc),
          plot_mt_genes(`fed4-1_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc`),
          plot_mt_genes(`fed4-2_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc`),
          plot_mt_genes(fed5_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc),
          plot_mt_genes(`fed6-1_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc`),
          plot_mt_genes(`fed6-2_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc`),
          # labels = "AUTO",
          # align = "hv",
          nrow = 3,
          ncol = 3)

plot_grid(plot_mt_genes(starv0_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc),
          plot_mt_genes(starv1_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc),
          plot_mt_genes(starv2_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc),
          plot_mt_genes(starv3_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc),
          plot_mt_genes(`starv4-2_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc`),
          plot_mt_genes(starv5_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc),
          plot_mt_genes(`starv6-1_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc`),
          plot_mt_genes(`starv6-2_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc`),
          # labels = "AUTO",
          # align = "hv",
          nrow = 3,
          ncol = 3)

# plot density plots of mt_genes percentage
plot_mt_genes_density <- function(df) {
  str <- deparse(substitute(df)) %>% str_split("_") %>% unlist() %>% .[1]
  p <-
    ggplot(df, aes(mt_genes/total)) +
    geom_density() +
    theme_classic() +
    theme(aspect.ratio = 1) +
    scale_x_continuous(limits = c(0, 1),
                       expand = expansion(mult = c(0, 0)),
                       breaks = seq(from = 0, to = 1, by = 0.1)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0))) +
    labs(title = str)
  
  return(p)
}

plot_grid(plot_mt_genes_density(fed0_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc),
          plot_mt_genes_density(fed1_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc),
          plot_mt_genes_density(fed2_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc),
          plot_mt_genes_density(fed3_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc),
          plot_mt_genes_density(`fed4-1_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc`),
          plot_mt_genes_density(`fed4-2_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc`),
          plot_mt_genes_density(fed5_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc),
          plot_mt_genes_density(`fed6-1_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc`),
          plot_mt_genes_density(`fed6-2_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc`),
          # labels = "AUTO",
          # align = "hv",
          nrow = 3,
          ncol = 3)

plot_grid(plot_mt_genes_density(starv0_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc),
          plot_mt_genes_density(starv1_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc),
          plot_mt_genes_density(starv2_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc),
          plot_mt_genes_density(starv3_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc),
          plot_mt_genes_density(`starv4-2_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc`),
          plot_mt_genes_density(starv5_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc),
          plot_mt_genes_density(`starv6-1_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc`),
          plot_mt_genes_density(`starv6-2_scKB_ws286_BestExtAgreedBy2reps_20230422_ByBc`),
          # labels = "AUTO",
          # align = "hv",
          nrow = 3,
          ncol = 3)

# combine samples and pick a cutoff (a stringent cut off to make sure the initial step of COPILOT QC picks only really bad barcodes)
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

p_mt_prop <- 
  ggplot(ws286_BestExt, aes(mt_prop)) +
  stat_ecdf() +
  theme_classic() +
  theme(aspect.ratio = 1) +
  scale_x_continuous(expand = expansion(mult = c(0, 0)),
                     breaks = seq(from = 0, to = 1, by = 0.05)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0)),
                     breaks = seq(from = 0, to = 1, by = 0.05)) +
  labs(title = "mt_prop all samples mapped to best ext")

p_mt_prop

quantile(ecdf(ws286_BestExt$mt_prop), 0.9) # 0.1666667
quantile(ecdf(ws286_BestExt$mt_prop), 0.95) # 0.25

# pick 0.2 mt_genes proportion as the really stringent cutoff

dev.off()


