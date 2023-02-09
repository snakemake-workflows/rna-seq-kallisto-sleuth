log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("sleuth")
library("tidyverse")

model <- snakemake@params[["model"]]

samples <- read_tsv(snakemake@input[["samples"]], na = "", col_names = TRUE) %>%
            # make everything except the index, sample name and path string a factor
            mutate_at(  vars(-sample, -path),
                        list(~factor(.))
                        )

if(!is.null(snakemake@params[["exclude"]])) {
    samples <- samples %>%
                filter( !sample %in% snakemake@params[["exclude"]] )
}

samples_out <- if(!is.null(model[["full"]])) {
    # retrieve the model formula
    formula <- as.formula(model[["full"]])
    # extract variables from the formula and unnest any nested variables
    variables <- labels(terms(formula)) %>%
                    strsplit('[:*]') %>%
                    unlist()
    # remove samples with an NA value in any of the columns
    # relevant for sleuth under the current model
    samples <- samples %>%
	        drop_na(c(sample, path, all_of(variables)))

    primary_variable <- model[["primary_variable"]]
    base_level <- model[["base_level"]]
    # TODO migrate this to tidyverse
    # Ensure that primary variable factors are sorted such that base_level comes first.
    # This is important for fold changes, effect sizes to have the expected sign.
    samples[, primary_variable] <- relevel(samples[, primary_variable, drop=TRUE], base_level)

    samples %>% select(c(sample, all_of(variables)))
} else {
    samples %>% select(-path)
}

# store design matrix
write_rds(samples_out, file = snakemake@output[["designmatrix"]])

# remove all columns which have only NA values
samples <- samples %>%
	    select_if(function(col) !all(is.na(col)))

t2g <- read_rds(snakemake@input[["transcript_info"]])

so <- sleuth_prep(  samples,
                    extra_bootstrap_summary = TRUE,
                    target_mapping = t2g,
                    aggregation_column = "ens_gene",
                    read_bootstrap_tpm = TRUE,
                    transform_fun_counts = function(x) log2(x + 0.5),
                    num_cores = snakemake@threads
                    )

if ("canonical" %in% colnames(so$target_mapping)) {
    # Sleuth converts all columns to characters. We don't want that for the canonical column.
    # Hence we fix it here.
    so$target_mapping$canonical <- as.logical(so$target_mapping$canonical)
}
custom_transcripts <- so$obs_raw %>%
                        # find transcripts not in the target_mapping
                        filter(!target_id %in% so$target_mapping$target_id) %>%
                        # make it a succinct list with no repititions
                        distinct(target_id) %>%
                        # pull it out into a vector
                        pull(target_id)

if(!length(custom_transcripts) == 0) {
    so$target_mapping <- so$target_mapping %>%
                        # add those custom transcripts as rows to the target mapping
                        add_row(ens_gene = NA, ext_gene = "Custom", target_id = custom_transcripts, canonical = NA)
}

if(!is.null(model[["full"]])) {
    so <- sleuth_fit(so, as.formula(model[["full"]]), 'full')
    so <- sleuth_fit(so, as.formula(model[["reduced"]]), 'reduced')
}

sleuth_save(so, snakemake@output[["sleuth_object"]])