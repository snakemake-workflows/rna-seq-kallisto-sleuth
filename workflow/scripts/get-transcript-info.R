log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("biomaRt")
library("tidyverse")

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
        version = snakemake@params[["version"]],
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
                            "external_gene_name",
                            "description"),
            mart = mart,
            useCache = FALSE
            ) %>%
        rename( target_id = ensembl_transcript_id,
                ens_gene = ensembl_gene_id,
                ext_gene = external_gene_name,
                gene_desc = description
                ) %>%
        mutate_at(
          vars(gene_desc),
          function(value) { str_trim(str_split(value, r"{\[}")[[1]][1]) } # remove trailing source annotation (e.g. [Source:HGNC Symbol;Acc:HGNC:5])
        )


write_rds(t2g, file = snakemake@output[[1]], compress = "gz")