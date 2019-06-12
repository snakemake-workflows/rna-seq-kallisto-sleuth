suppressMessages({
  library("sleuth")
  library("biomaRt")
})

model <- snakemake@params[["model"]]

mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = paste(snakemake@params[["species"]], "_gene_ensembl", sep = ""),
  host = 'ensembl.org')
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
    "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
  ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

samples <- read.table(snakemake@input[["samples"]], colClasses = "factor", sep = "\t", na.strings = "", header = TRUE)
samples[, "path"] <- as.character(samples[, "path"])

if(!is.null(snakemake@params[["exclude"]])) {
    samples <- dplyr::filter(samples, !samples$sample %in% snakemake@params[["exclude"]] )
}

if(!is.null(model)) {
    formula <- as.formula(snakemake@params[["model"]])
    variables <- colnames(attr(terms(formula), "factors"))
    cols <- c("sample", "path", "condition", unlist(strsplit(variables, '[:*]')) )
    # filter to all samples where the given variables are not NA
    samples <- samples[complete.cases(samples[, cols]), cols]
}

so <- sleuth_prep(samples, extra_bootstrap_summary = TRUE, target_mapping = t2g, aggregation_column = "ens_gene")

# add custom transcripts to target mapping
custom_transcripts <- unique(so$obs_raw[(!so$obs_raw$target_id %in% so$target_mapping$target_id), "target_id"])
so$target_mapping <- rbind(so$target_mapping, data.frame(ens_gene = NA, ext_gene = "Custom", target_id = custom_transcripts))

sleuth_save(so, snakemake@output[[1]])
