log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("sleuth")
library("biomaRt")
library("tidyverse")

model <- snakemake@params[["model"]]

mart <- biomaRt::useMart(
            biomart = "ENSEMBL_MART_ENSEMBL",
            dataset = str_c(snakemake@params[["species"]], "_gene_ensembl"),
            host = 'ensembl.org'
            )
t2g <- biomaRt::getBM(
            attributes = c( "ensembl_transcript_id",
                            "ensembl_gene_id",
                            "external_gene_name"),
            mart = mart
            ) %>%
        rename( target_id = ensembl_transcript_id,
                ens_gene = ensembl_gene_id,
                ext_gene = external_gene_name
                )

samples <- read_tsv(snakemake@input[["samples"]], na = "", col_names = TRUE) %>%
            # make everything except the index, sample name and path string a factor
            mutate_at(  vars(-X1, -sample, -path),
                        list(~factor(.))
                        )

if(!is.null(snakemake@params[["exclude"]])) {
    samples <- samples %>%
                filter( !sample %in% snakemake@params[["exclude"]] )
}

if(!is.null(model)) {
    # retrieve the model formula
    formula <- as.formula(snakemake@params[["model"]])
    # extract variables from the formula and unnest any nested variables
    variables <- labels(terms(formula)) %>%
                    strsplit('[:*]') %>%
                    unlist()
    # select the columns required by sleuth and filter to all samples where
    # none of the given variables are NA
    samples <- samples %>%
                select(sample, path, condition, variables) %>%
                drop_na()
}

so <- sleuth_prep(  samples,
                    extra_bootstrap_summary = TRUE,
                    target_mapping = t2g,
                    aggregation_column = "ens_gene",
                    read_bootstrap_tpm = TRUE,
                    transform_fun_counts = function(x) log2(x + 0.5),
                    num_cores = snakemake@threads
                    )

custom_transcripts <- so$obs_raw %>%
                        # find transcripts not in the target_mapping
                        filter(!target_id %in% so$target_mapping$target_id) %>%
                        # make it a succinct list with no repititions
                        distinct(target_id) %>%
                        # pull it out into a vector
                        pull(target_id)
so$target_mapping <- so$target_mapping %>%
                        # add those custom transcripts as rows to the target mapping
                        add_row(ens_gene = NA, ext_gene = "Custom", target_id = custom_transcripts)

sleuth_save(so, snakemake@output[[1]])
