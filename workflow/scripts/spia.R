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
        "Name", "Status", "Combined FDR",
        "total perturbation accumulation", "number of genes on the pathway",
        "number of DE genes per pathway", "p-value for at least NDE genes",
        "Combined Bonferroni p-values",
        "p-value to observe a total accumulation", "Combined p-value", "Ids"
    )
    res <- data.frame(matrix(ncol = 11, nrow = 0, dimnames = list(NULL, cols)))
    # create empty perturbation plots
    pdf(file = snakemake@output[["plots"]])
    write_tsv(res, snakemake@output[["table"]])
    write_tsv(res, snakemake@output[["table_activated"]])
    write_tsv(res, snakemake@output[["table_inhibited"]])
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
                "Combined Bonferroni p-values" = "pGFWER",
                "Combined FDR" = "pGFdr",
                "total perturbation accumulation" = "tA",
                "number of genes on the pathway" = "pSize",
                "number of DE genes per pathway" = "NDE",
                "Combined p-value" = "pG",
                "p-value to observe a total accumulation" = "pPERT",
                "p-value for at least NDE genes" = "pNDE"
            )
        write_tsv(res_reorder, snakemake@output[["table"]])
        sort_activated <- res_reorder[res_reorder$Status == "Activated", ]
        sort_inhibited <- res_reorder[res_reorder$Status == "Inhibited", ]
        write_tsv(sort_activated, snakemake@output[["table_activated"]])
        write_tsv(sort_inhibited, snakemake@output[["table_inhibited"]])
    } else {
        columns <- c(
            "Name", "Status", "Combined FDR", "total perturbation accumulation",
            "number of genes on the pathway", "number of DE genes per pathway",
            "p-value for at least NDE genes", "Combined Bonferroni p-values",
            "p-value to observe a total accumulation", "Combined p-value", "Ids"
        )
        emtpy_data_frame <- data.frame(matrix(nrow = 0, ncol = length(columns)))
        colnames(emtpy_data_frame) <- columns
        write_tsv(emtpy_data_frame, snakemake@output[["table"]])
        write_tsv(emtpy_data_frame, snakemake@output[["table_activated"]])
        write_tsv(emtpy_data_frame, snakemake@output[["table_inhibited"]])
    }
}
