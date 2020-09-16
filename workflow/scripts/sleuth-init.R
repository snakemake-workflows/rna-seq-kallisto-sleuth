log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("sleuth")
library("biomaRt")
library("tidyverse")

model <- snakemake@params[["model"]]

# this variable holds a mirror name until
# useEnsembl succeeds ("www" is last, because 
# of very frequent "Internal Server Error"s)
mart <- "useast"
rounds <- 0
while ( class(mart)[[1]] != "Mart" ) {
  mart <- tryCatch(
    {
      # done here, because error function does not
      # modify outer scope variables, I tried
      if (mart == "www") rounds <- rounds + 1
      # equivalent to useMart, but you can choose
      # the mirror instead of specifying a host
      biomaRt::useEnsembl(
        biomart = "ENSEMBL_MART_ENSEMBL",
        dataset = str_c(snakemake@params[["species"]], "_gene_ensembl"),
        mirror = mart
      )
    },
    error = function(e) {
      # change or make configurable if you want more or
      # less rounds of tries of all the mirrors
      if (rounds >= 3) {
        stop(
          str_c(
            "Have tried all 4 available Ensembl biomaRt mirrors ",
            rounds,
            " times. You might have a connection problem, or no mirror is responsive."
          )
        )
      }
      # hop to next mirror
      mart <- switch(mart,
                     useast = "uswest",
                     uswest = "asia",
                     asia = "www",
                     www = {
                       # wait before starting another round through the mirrors,
                       # hoping that intermittent problems disappear
                       Sys.sleep(30)
                       "useast"
                     }
              )
    }
  )
}

t2g <- biomaRt::getBM(
            attributes = c( "ensembl_transcript_id",
                            "ensembl_gene_id",
                            "external_gene_name"),
            mart = mart,
            useCache = FALSE
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
    formula <- as.formula(model)
    # extract variables from the formula and unnest any nested variables
    variables <- labels(terms(formula)) %>%
                    strsplit('[:*]') %>%
                    unlist()
    # remove samples with an NA value in any of the columns
    # relevant for sleuth under the current model
    samples <- samples %>%
                drop_na(
                    c( sample, path, all_of(variables) )
                )
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

if(!length(custom_transcripts) == 0) {
    so$target_mapping <- so$target_mapping %>%
                        # add those custom transcripts as rows to the target mapping
                        add_row(ens_gene = NA, ext_gene = "Custom", target_id = custom_transcripts)
}

sleuth_save(so, snakemake@output[[1]])
