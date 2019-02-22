suppressMessages({
  library("sleuth")
})

so <- sleuth_load(snakemake@input[[1]])

so <- sleuth_fit(so, as.formula(snakemake@params[["model"]]), 'full')
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')

all <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE, pval_aggregate = FALSE)
all <- dplyr::arrange(all, pval)

all_aggregated <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE, pval_aggregate = TRUE)
all_aggregated <- dplyr::arrange(all_aggregated, pval)

write.table(all, file = snakemake@output[[1]], quote=FALSE, sep='\t', row.names = FALSE)
write.table(all_aggregated, snakemake@output[[2]], quote=FALSE, sep='\t', row.names = FALSE)
