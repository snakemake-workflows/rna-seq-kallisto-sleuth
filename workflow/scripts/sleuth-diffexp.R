log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("sleuth")
library("tidyverse")
library("fs")

model <- snakemake@params[["model"]]

dir_create(file.path(snakemake@output[["volcano_plots"]]))

sleuth_object <- sleuth_load(snakemake@input[[1]])

sleuth_object <- sleuth_fit(sleuth_object, as.formula(model[["full"]]), 'full')
sleuth_object <- sleuth_fit(sleuth_object, as.formula(model[["reduced"]]), 'reduced')
sleuth_object <- sleuth_lrt(sleuth_object, "reduced", "full")

write_results <- function(so, mode, output, output_all) {
    so$pval_aggregate <- FALSE
    if(mode == "aggregate") {
      # workaround the following bug-request:
      # https://github.com/pachterlab/sleuth/pull/240
      # TODO renaming can be removed when fixed
      g_col <- so$gene_column
      so$gene_column <- NULL
      so$pval_aggregate <- TRUE
      so$gene_column <- g_col
    }

    print("Performing likelihood ratio test")
    all <- sleuth_results(so, "reduced:full", "lrt", show_all = TRUE, pval_aggregate = so$pval_aggregate) %>%
            arrange(pval)

    covariates <- colnames(design_matrix(so, "full"))
    covariates <- covariates[covariates != "(Intercept)"]

    # iterate over all covariates and perform wald test in order to obtain beta estimates
    if(!so$pval_aggregate) {
        for(covariate in covariates) { 
            print(str_c("Performing wald test for ", covariate))
            so <- sleuth_wt(so, covariate, "full")

            print(str_c("Performing volcano plot for ", covariate))
            path_output <- str_c(snakemake@output[["volcano_plots"]], "/", snakemake@wildcards[["model"]],".volcano-plot.", covariate, ".pdf")
            
            volcano <- plot_volcano(so, covariate, "wt", "full",
                                         sig_level = 0.1, point_alpha = 0.2, sig_color = "red",
                                         highlight = NULL)
 
            ggsave(path_output, plot = volcano, width = 14)

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
write_results(sleuth_object, "transcripts", snakemake@output[["transcripts"]], snakemake@output[["transcripts_rds"]])
write_results(sleuth_object, "aggregate", snakemake@output[["genes_aggregated"]], snakemake@output[["genes_aggregated_rds"]])
write_results(sleuth_object, "mostsignificant", snakemake@output[["genes_mostsigtrans"]], snakemake@output[["genes_mostsigtrans_rds"]])
