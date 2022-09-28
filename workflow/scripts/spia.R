log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("SPIA")
library("graphite")

# provides library("tidyverse") and functions load_bioconductor_package() and
# get_prefix_col(), the latter requires snakemake@output[["samples"]] and
# snakemake@params[["covariate"]]
source(snakemake@params[["common_src"]])

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

if(nrow(sig_genes) == 0) {
    cols <- c("Name", "pSize", "NDE", "pNDE", "tA", "pPERT", "pG", "pGFdr", "pGFWER", "Status")
    res <- data.frame(matrix(ncol=10, nrow=0, dimnames=list(NULL, cols)))
    # create empty perturbation plots
    pdf(file = snakemake@output[["plots"]])
    dev.off()
} else {

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
    pathway_names <- db[res$Name]
    pathway_names <- db[res$Name]
    path_ids <- as.matrix(lapply(pathway_names@entries, slot, "id"))
    path_ids_data_frame <-
        data.frame(Ids = matrix(unlist(path_ids),
            nrow = length(path_ids), byrow = TRUE))
    final_res <- cbind(res,
        Ids = path_ids_data_frame$Ids)
    res_reorder <- dplyr::select(final_res, Status, Name,
        pGFdr, tA, pSize, NDE, pNDE, pGFWER, pPERT, pG, Ids)
    names(res_reorder)[names(res_reorder) == "pGFWER"] <- "Bonferroni adjusted global p-values"
    names(res_reorder)[names(res_reorder) == "pGFdr"] <- "False Discovery Rate global p-values"
    names(res_reorder)[names(res_reorder) == "tA"] <- "observed total perturbation accumulation"
    names(res_reorder)[names(res_reorder) == "pSize"] <- "number of genes on the pathway"
    names(res_reorder)[names(res_reorder) == "NDE"] <- "number of DE genes per pathway"
    names(res_reorder)[names(res_reorder) == "pG"] <- "p-value obtained by combining pNDE and pPERT"
    names(res_reorder)[names(res_reorder) == "pPERT"] <- "probability to observe a total accumulation"
    names(res_reorder)[names(res_reorder) == "pNDE"] <- "probability to observe at least NDE genes on the pathway"
    write_tsv(res_reorder, snakemake@output[["table"]])
}