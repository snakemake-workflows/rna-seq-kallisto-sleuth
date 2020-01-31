log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("SPIA")
library("graphite")

# provides library("tidyverse") and functions load_bioconductor_package() and
# get_prefix_col(), the latter requires snakemake@output[["samples"]] and
# snakemake@params[["covariate"]]
source( file.path(snakemake@scriptdir, 'common.R') )

pkg <- snakemake@params[["bioc_pkg"]]
load_bioconductor_package(snakemake@input[["species_anno"]], pkg)



pw_db <- snakemake@params[["pathway_db"]]

db <- pathways(snakemake@params[["species"]], pw_db)
db <- convertIdentifiers(db, "ENSEMBL")

options(Ncpus = snakemake@threads)

diffexp <- read_tsv(snakemake@input[["diffexp"]]) %>%
            drop_na(ens_gene) %>%
            mutate(ens_gene = str_c("ENSEMBL:", ens_gene))
universe <- diffexp %>% pull(var = ens_gene)
sig_genes <- diffexp %>% filter(qval <= 0.05)

# get logFC equivalent (the sum of beta scores of covariates of interest)

beta_col <- get_prefix_col("b", colnames(sig_genes))

beta <- sig_genes %>%
            dplyr::select(ens_gene, !!beta_col) %>%
            deframe()

t <- tempdir(check=TRUE)
olddir <- getwd()
setwd(t)
prepareSPIA(db, pw_db)
res <- runSPIA(de = beta, all = universe, pw_db, plots = TRUE, verbose = TRUE)
setwd(olddir)

file.copy(file.path(t, "SPIAPerturbationPlots.pdf"), snakemake@output[["plots"]])
write_tsv(res, snakemake@output[["table"]])
