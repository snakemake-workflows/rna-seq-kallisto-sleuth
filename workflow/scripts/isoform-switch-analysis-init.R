log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("tidyverse")
library("sleuth")
library("IsoformSwitchAnalyzeR")

model <- snakemake@params[["model"]]

m <- readRDS(snakemake@input[["designmatrix"]])
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
    reduceToSwitchingGenes = TRUE,
    dIFcutoff = snakemake@params[["min_effect_size"]],
    alpha = snakemake@params[["fdr"]]
)

extractSequence(
    results,
    pathToOutput = snakemake@params[["seq_dir"]],
    addToSwitchAnalyzeRlist = TRUE
)

saveRDS(results, file = snakemake@output[[1]])
