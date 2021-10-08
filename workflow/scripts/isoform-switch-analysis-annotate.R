log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("tidyverse")
library("IsoformSwitchAnalyzeR")

results <- read_rds(snakemake@input[["rds"]])

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

results <- analyzeAlternativeSplicing(
    results, 
    quiet = FALSE, 
    onlySwitchingGenes = FALSE,
)

dir.create(snakemake@output[["plots_with"]])
dir.create(snakemake@output[["plots_without"]])

if(nrow(results$isoformFeatures) > 0) {
    results <- analyzeSwitchConsequences(
        results, 
        consequencesToAnalyze = c(
            'intron_retention',
            'coding_potential',
            # 'ORF_seq_similarity', TODO this is only needed for assembly, reactivate then
            'NMD_status',
            'domains_identified'
        ),
        onlySigIsoforms = FALSE,
        removeNonConseqSwitches = FALSE,
        quiet = TRUE, # set to TRUE to circumvent a bug leading to a stop() if there are no significant switches
        alpha = snakemake@params[["fdr"]],
        dIFcutoff = snakemake@params[["min_effect_size"]],
    )

    switchPlotTopSwitches(
        switchAnalyzeRlist = results,
        n = Inf,
        filterForConsequences = FALSE,
        splitComparison = FALSE,
        splitFunctionalConsequences = TRUE,
        pathToOutput = snakemake@params[["plotdir"]],
        alpha = snakemake@params[["fdr"]],
        dIFcutoff = snakemake@params[["min_effect_size"]],
        onlySigIsoforms = FALSE,
    )
}

significant <- extractTopSwitches(
    results,
    filterForConsequences = FALSE,
    extractGenes = TRUE,
    n = Inf,
    sortByQvals = TRUE,
    alpha = snakemake@params[["fdr"]],
    dIFcutoff = snakemake@params[["min_effect_size"]],
)

write_tsv(significant, snakemake@output[["table"]])
