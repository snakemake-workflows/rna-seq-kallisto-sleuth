suppressMessages({
  library("sleuth")
})

model <- snakemake@params[["model"]]

so <- sleuth_load(snakemake@input[[1]])

so <- sleuth_fit(so, as.formula(model[["full"]]), 'full')
so <- sleuth_fit(so, as.formula(model[["reduced"]]), 'reduced')
so <- sleuth_lrt(so, "reduced", "full")
so <- sleuth_wt(so, beta_covariate, "full")

write_results <- function(aggregate, output) {
    all <- sleuth_results(so, "reduced:full", "lrt", show_all = FALSE, pval_aggregate = aggregate)
    all <- dplyr::arrange(all, pval)
    for(covariate in colnames(design_matrix(so, "full"))) {	    
        all_wald <- sleuth_results(so, covariate, "wt", show_all = FALSE, pval_aggregate = aggregate)
        all_wald <- all_wald[, c("target_id", "b", "se_b")]
	colnames(all_wald) <- c("target_id", paste("b", covariate, sep = "_"), paste("b", covariate, "se", sep = "_"))
        all <- merge(all, all_wald, on = "target_id", sort = FALSE)
    }
    write.table(all, file = snakemake@output[["transcripts"]], quote=FALSE, sep='\t', row.names = FALSE)
}

write_results(FALSE, snakemake@output[["transcripts"]])
write_results(TRUE, snakemake@output[["genes"]])

