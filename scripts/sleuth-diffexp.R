suppressMessages({
  library("sleuth")
})

model <- snakemake@params[["model"]]

so <- sleuth_load(snakemake@input[[1]])

so <- sleuth_fit(so, as.formula(model[["full"]]), 'full')
so <- sleuth_fit(so, as.formula(model[["reduced"]]), 'reduced')
so <- sleuth_lrt(so, "reduced", "full")

write_results <- function(mode, output, output_all) {
    aggregate <- FALSE
    if(mode == "aggregate") {
        aggregate <- TRUE
    }
    print("Performing likelihood ratio test")
    all <- sleuth_results(so, "reduced:full", "lrt", show_all = FALSE, pval_aggregate = aggregate)
    all <- dplyr::arrange(all, pval)

    covariates <- colnames(design_matrix(so, "full"))
    covariates <- covariates[covariates != "(Intercept)"]

    # iterate over all covariates and perform wald test in order to obtain beta estimates
    if(!aggregate) {
        for(covariate in covariates) { 
            print(paste("Performing wald test for", covariate))
            so <- sleuth_wt(so, covariate, "full")

            all_wald <- sleuth_results(so, covariate, "wt", show_all = TRUE, pval_aggregate = FALSE)
            all_wald <- all_wald[, c("target_id", "b", "se_b")]
	    beta <- paste("b", covariate, sep = "_")
            colnames(all_wald) <- c("target_id", beta, paste(beta, "se", sep = "_"))

            all <- merge(all, all_wald, on = "target_id", sort = FALSE)
	    # calculate a signed version of the pi-value from https://dx.doi.org/10.1093/bioinformatics/btr671
	    # e.g. useful for GSEA ranking
	    all[, paste("signed_pi_value", covariate, sep = "_")] <- -log10(all$pval) * all[, beta]
        }
    }

    if(mode == "mostsignificant") {
        # for each gene, find most significant (this is equivalent to sleuth_gene_table, but with bug fixes)
        gene_table <- dplyr::arrange_(gene_table, "ext_gene", ~qval)
        gene_table <- dplyr::group_by_(gene_table, "ext_gene")
	gene_table <- dplyr::summarise_(gene_table, target_id = ~target_id[1], pval = ~min(pval, na.rm = TRUE), qval = ~min(qval, na.rm = TRUE),  num_transcripts = ~n(), all_target_ids = ~paste0(target_id[1:length(target_id)], collapse = ','))
        # only keep those transcripts
	all <- all[gene_table$target_id, ]

    }

    write.table(all, file = output, quote=FALSE, sep='\t', row.names = FALSE)
    saveRDS(all, file = output_all, compress = FALSE)
}

write_results("transcripts", snakemake@output[["transcripts"]], snakemake@output[["transcripts_rds"]])
write_results("aggregate", snakemake@output[["genes_aggregated"]], snakemake@output[["genes_aggregated_rds"]])
write_results("mostsigificant", snakemake@output[["genes_mostsigtrans"]], snakemake@output[["genes_mostsigtrans_rds"]])
