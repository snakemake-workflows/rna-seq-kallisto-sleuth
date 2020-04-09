log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("tidyverse")
library("IsoformSwitchAnalyzeR")

results <- readRDS(snakemake@input[["rds"]])

results <- analyzePFAM(
    switchAnalyzeRlist = results, 
    pathToPFAMresultFile = snakemake@input[["pfam"]],
    quiet = TRUE
)

results <- analyzeCPAT(
    switchAnalyzeRlist = results,
    pathToCPATresultFile = snakemake@input[["cpat"]],
    codingCutoff = snakemake@params[["coding_cutoff"]],
    removeNoncodinORFs = snakemake@params[["remove_noncoding_orfs"]],
    quiet = TRUE
)

results <- analyzeAlternativeSplicing(results, quiet = FALSE)

results <- analyzeSwitchConsequences(
    results, 
    consequencesToAnalyze = c(
        'intron_retention',
	'coding_potential',
	# 'ORF_seq_similarity', TODO this is only needed for assembly, reactivate then
	'NMD_status',
	'domains_identified'
    ),
    quiet = FALSE
)

switchPlotTopSwitches(
    switchAnalyzeRlist = results,
    n = Inf,
    filterForConsequences = FALSE,
    splitComparison = FALSE,
    splitFunctionalConsequences = TRUE,
    pathToOutput = snakemake@params[["plotdir"]],
)

significant <- extractTopSwitches(
    results,
    filterForConsequences = FALSE,
    extractGenes = TRUE,
    n = Inf,
    sortByQvals = TRUE,
)

write_tsv(significant, snakemake@output[["table"]])
