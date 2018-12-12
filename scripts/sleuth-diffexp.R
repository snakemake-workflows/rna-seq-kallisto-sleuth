suppressMessages({
  library("sleuth")
})

so <- sleuth_load(snakemake@input[[1]])

so <- sleuth_fit(so, as.formula(snakemake@params[["model"]]), 'full')
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')

all <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
all <- dplyr::arrange(all, pval)

write.table(all, file = snakemake@output[[1]], quote=FALSE, sep='\t', row.names = FALSE)
