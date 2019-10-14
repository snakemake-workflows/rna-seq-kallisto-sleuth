log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("sleuth")
library("tidyverse")

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
    all <- sleuth_results(so, "reduced:full", "lrt", show_all = TRUE, pval_aggregate = aggregate) %>%
            arrange(pval)

    covariates <- colnames(design_matrix(so, "full"))
    covariates <- covariates[covariates != "(Intercept)"]

    # iterate over all covariates and perform wald test in order to obtain beta estimates
    if(!aggregate) {
        for(covariate in covariates) { 
            print(str_c("Performing wald test for ", covariate))
            so <- sleuth_wt(so, covariate, "full")

	      beta_col_name <- str_c("b", covariate, sep = "_")
            beta_se_col_name <- str_c(beta_col_name, "se", sep = "_")
            all_wald <- sleuth_results(so, covariate, "wt", show_all = TRUE, pval_aggregate = FALSE) %>%
                        select( target_id = target_id,
                                !!beta_col_name := b,
                                !!beta_se_col_name := se_b)
            signed_pi_col_name <- str_c("signed_pi_value", covariate, sep = "_")
            all <- inner_join(all, all_wald, by = "target_id") %>%
	              # calculate a signed version of the pi-value from:
                    # https://dx.doi.org/10.1093/bioinformatics/btr671
	              # e.g. useful for GSEA ranking
	              mutate( !!signed_pi_col_name := -log10(pval) * !!sym(beta_col_name) )
        }
    }

    if(mode == "mostsignificant") {
        # for each gene, select the most significant transcript (this is equivalent to sleuth_gene_table, but with bug fixes)
	  all <- all %>%
                drop_na(ens_gene) %>%
                group_by(ens_gene) %>%
                filter( qval == min(qval, na.rm = TRUE) ) %>%
                # ties in qval (e.g. min(qval) == 1) can lead to multiple entries per gene
                filter( pval == min(pval, na.rm = TRUE) ) %>%
                # for min(qval) == 1, then usually also min(pval) == 1, and
                # we need something else to select a unique entry per gene
                filter( mean_obs == max(mean_obs, na.rm = TRUE) ) %>%
                # sometimes we still get two transcript variants with exactly
                # the same values, i.e. they have exactly the same sequence
                # but (slightly) different annotations -- then we retain a string
                # with a comma-separated list of them
                mutate(target_id = str_c(target_id, collapse=",")) %>%
                distinct() %>%
                # useful sort for scrolling through output by increasing q-values
                arrange(qval)
    }

    write_tsv(all, path = output, quote_escape = "none")
    write_rds(all, path = output_all, compress = "none")
}

write_results("transcripts", snakemake@output[["transcripts"]], snakemake@output[["transcripts_rds"]])
write_results("aggregate", snakemake@output[["genes_aggregated"]], snakemake@output[["genes_aggregated_rds"]])
write_results("mostsignificant", snakemake@output[["genes_mostsigtrans"]], snakemake@output[["genes_mostsigtrans_rds"]])
