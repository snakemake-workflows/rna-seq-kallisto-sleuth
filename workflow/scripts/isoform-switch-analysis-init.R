log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("tidyverse")
library("IsoformSwitchAnalyzeR")

model <- snakemake@params[["model"]]

m <- read_rds(snakemake@input[["designmatrix"]])
is_prefix_col <- startsWith(colnames(m), model[["primary_variable"]])
colnames(m) <- replace(colnames(m), is_prefix_col, "condition")
colnames(m) <- replace(colnames(m), colnames(m) == "sample", "sampleID")
m <- m %>%
      dplyr::select("sampleID", "condition", everything())

quant <- importIsoformExpression(
   sampleVector = set_names(paste(snakemake@input[["kallisto"]], "abundance.tsv", sep="/"), snakemake@params[["samples"]]),
   addIsofomIdAsColumn = TRUE
)

candidates <- importRdata(
    isoformCountMatrix = quant$counts,
    isoformRepExpression = quant$abundance,
    designMatrix = m,
    isoformExonAnnoation = snakemake@input[["gtf"]],
    isoformNtFasta = snakemake@input[["fasta"]],
    showProgress = FALSE
)

# TODO make cutoffs configurable
filtered_candidates <- preFilter(
    switchAnalyzeRlist = candidates,
    geneExpressionCutoff = 1,
    isoformExpressionCutoff = 0,
    removeSingleIsoformGenes = TRUE
)

results <- isoformSwitchTestDEXSeq(
    switchAnalyzeRlist = filtered_candidates,
    reduceToSwitchingGenes = FALSE,
)

# get significant genes (have to do this manually, because IsoformSwitchAnalyzeR exits with an error
# in case of no significant genes).
keep <- unique(
    results$isoformFeatures$gene_ref[which(
        results$isoformFeatures$isoform_switch_q_value <
            snakemake@params[["fdr"]] &
            abs(results$isoformFeatures$dIF) > snakemake@params[["min_effect_size"]]
    )]
)

# if(length(keep) == 0) {
#     # No significant genes left, just keep the first to keep the analysis going.
#     keep <- results$isoformFeatures$gene_ref[1]
# }

# filter to significant genes
results <- subsetSwitchAnalyzeRlist(results, results$isoformFeatures$gene_ref %in% keep)

extractSequence(
    results,
    pathToOutput = snakemake@params[["seq_dir"]],
    addToSwitchAnalyzeRlist = TRUE,
    onlySwitchingGenes = FALSE,
)

write_rds(results, file = snakemake@output[[1]], compress = "gz")
