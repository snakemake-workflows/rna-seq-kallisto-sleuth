log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library("SPIA")
library("graphite")
library(snakemake@params[["bioc_species_pkg"]], character.only = TRUE)

# provides library("tidyverse") and functions load_bioconductor_package() and
# get_prefix_col(), the latter requires snakemake@output[["samples"]] and
# snakemake@params[["covariate"]]
source(snakemake@params[["common_src"]])

pw_db <- snakemake@params[["pathway_db"]]
db <- readRDS(snakemake@input[["spia_db"]])

options(Ncpus = snakemake@threads)

diffexp <- read_tsv(snakemake@input[["diffexp"]]) %>%
    drop_na(ens_gene) %>%
    mutate(ens_gene = str_c("ENSEMBL:", ens_gene))
universe <- diffexp %>% pull(var = ens_gene)
sig_genes <- diffexp %>% filter(qval <= 0.05)

if (nrow(sig_genes) == 0) {
    cols <- c(
        "Name", "pSize", "NDE", "pNDE", "tA",
        "pPERT", "pG", "pGFdr", "pGFWER", "Status"
    )
    res <- data.frame(matrix(ncol = 10, nrow = 0, dimnames = list(NULL, cols)))
    # create empty perturbation plots
    pdf(file = snakemake@output[["plots"]])
    dev.off()
} else {
    # get logFC equivalent (the sum of beta scores of covariates of interest)

    beta_col <- get_prefix_col("b", colnames(sig_genes))

    beta <- sig_genes %>%
        dplyr::select(ens_gene, !!beta_col) %>%
        deframe()

    t <- tempdir(check = TRUE)
    olddir <- getwd()
    setwd(t)
    prepareSPIA(db, pw_db)
    res <- runSPIA(
        de = beta, all = universe, pw_db,
        plots = TRUE, verbose = TRUE
    )
    setwd(olddir)

    file.copy(
        file.path(t, "SPIAPerturbationPlots.pdf"),
        snakemake@output[["plots"]]
    )
    pathway_names <- db[res$Name]
    pathway_names <- db[res$Name]
    path_ids <- as.matrix(lapply(pathway_names@entries, slot, "id"))
    if (length(path_ids) > 0) {
        path_ids_data_frame <-
            data.frame(Ids = matrix(unlist(path_ids),
                nrow = length(path_ids), byrow = TRUE
            ))
        final_res <- cbind(res,
            Ids = path_ids_data_frame$Ids
        )
        res_reorder <- dplyr::select(
            final_res, Name, Status,
            pGFdr, tA, pSize, NDE, pNDE, pGFWER, pPERT, pG, Ids
        )
        res_reorder <- res_reorder %>%
            rename(
                "pGFdr" = "Combined Bonferroni p-values",
                "tA" = "total perturbation accumulation",
                "pSize" = "number of genes on the pathway",
                "NDE" = "number of DE genes per pathway",
                "pG" = "Combined p-value",
                "pPERT" = "p-value to observe a total accumulation",
                "pNDE" = "p-value for at least NDE genes"
            )
        write_tsv(res_reorder, snakemake@output[["table"]])
    } else {
        columns <- c(
            "Combined Bonferroni p-values", "Combined FDR",
            "total perturbation accumulation", "number of genes on the pathway",
            "Combined p-value no", "p-value to observe a total accumulation",
            "p-value for at least NDE genes"
        )
        emtpy_data_frame <- data.frame(matrix(nrow = 0, ncol = length(columns)))
        colnames(emtpy_data_frame) <- columns
        write_tsv(emtpy_data_frame, snakemake@output[["table"]])
    }
}
