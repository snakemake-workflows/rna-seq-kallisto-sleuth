suppressMessages({
  library("sleuth")
})

print(snakemake@params[["samples"]])

samples <- read.table(snakemake@input[["samples"]])
samples[, "path"] <- as.character(samples[, "path"])
print(samples)

so <- sleuth_prep(samples, extra_bootstrap_summary = TRUE)

sleuth_save(so, snakemake@output[[1]])
