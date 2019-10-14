log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("sleuth")
library("ggplot2")
library("tidyr")

diffexp <- sleuth_load(snakemake@input[["diffexp_rds"]]) %>% drop_na(pval)

ggplot(diffexp) + geom_histogram(aes(pval), bins = 100)
ggsave(file = snakemake@output[[1]], width = 14)
