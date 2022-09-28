log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(tidyverse)
library(ggplot2)


pathway_file <- read.csv(snakemake@input[["pathway_file"]],
    sep = "\t", header = TRUE)
fil_pathway <- pathway_file %>% filter(False.Discovery.Rate.global.p.values <= 0.05)
pdf(file = snakemake@output[[1]], height = 10, width = 10)
ggplot(fil_pathway, aes(reorder(Name, observed.total.perturbation.accumulation),
    observed.total.perturbation.accumulation)) +
    geom_col(aes(fill = Status)) + coord_flip() +
    labs(x = "Pathway",
        y = "the observed total perturbation accumulation in the pathway",
            title = "pathways tA from SPIA of pGfdr<05") +
    theme_minimal()
dev.off()
