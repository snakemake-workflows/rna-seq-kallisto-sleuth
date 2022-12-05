log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("SPIA")
library("graphite")
library(snakemake@params[["bioc_species_pkg"]], character.only = TRUE)

# provides library("tidyverse") and functions load_bioconductor_package() and
# get_prefix_col(), the latter requires snakemake@output[["samples"]] and
# snakemake@params[["covariate"]]
source(snakemake@params[["common_src"]])

pw_db <- snakemake@params[["pathway_db"]]

db <- pathways(snakemake@params[["species"]], pw_db)
db <- convertIdentifiers(db, "ENSEMBL")

saveRDS(db, snakemake@output[[1]])