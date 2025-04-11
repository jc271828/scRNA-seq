library(tidyverse)
library(magrittr)
library(readxl)
library(writexl)

bulk <- read_xlsx("JC_N2starv_vs_N2fed_RNASeq_EtOH_16hPostHatch.xlsx")
sc <- read_xlsx("fed_starv_compr_SCT.xlsx")

sc$gene_id %<>% str_remove("\\.")
sc %<>% mutate(log2FC_sc = log2(starv_SCT_avg/fed_SCT_avg))
colnames(bulk)[12] <- "log2FC_bulk"

compr <- merge(bulk, sc, by = "gene_id", all = F)

compr <- compr[is.finite(compr$log2FC_sc), ]

ggplot() +
  stat_function(color = "red", linewidth = 1.5, fun = function(x) x) +
  geom_vline(xintercept = 0, color = "red", linewidth = 1.5) +
  geom_point(data = filter(compr, abs(log2FC_sc) <= 2),
             color = "gray",
             mapping = aes(x = log2FC_sc,
                           y = log2FC_bulk)) +
  geom_point(data = filter(compr, abs(log2FC_sc) > 2),
             mapping = aes(x = log2FC_sc,
                           y = log2FC_bulk)) +
  theme_classic() +
  theme(aspect.ratio = 0.5) +
  # coord_cartesian(ylim = c(0, 10)) +
  scale_y_continuous(limits = c(0, 10), expand = expansion(mult = c(0, 0))) +
  scale_x_continuous(limits = c(-10, 10), expand = expansion(mult = c(0, 0))) +
  labs(title = "starvation response genes log2FC changes in sc vs bulk RNASeq\nscRNASeq data is non-proportional pseudo bulk\ngene list from JC's N2 starv vs fed RNASeq @ 16h post hatch in 0.1% EtOH\ngray: log2FC_sc absolute value <= 2",
       x = "log2FC in scRNASeq (non-proportional pseudo bulk)",
       y = "log2FC in bulk RNASeq")

ggplot(data = compr) +
  stat_ecdf(mapping = aes(x = log2FC_sc), linewidth = 1) +
  theme_classic() +
  theme(aspect.ratio = 1)


