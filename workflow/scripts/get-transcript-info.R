 log <- file(snakemake@log[[1]], open="wt")
 sink(log)
 sink(log, type="message")

library("biomaRt")
library("tidyverse")
library("dplyr")

# this variable holds a mirror name until
# useEnsembl succeeds ("www" is last, because
# of very frequent "Internal Server Error"s)
mart <- "useast"
rounds <- 0
while (class(mart)[[1]] != "Mart") {
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
three_prime_activated <- snakemake@params[["three_prime_activated"]]

attributes <- c("ensembl_transcript_id",
                "ensembl_gene_id",
                "external_gene_name",
                "description")
has_canonical <-
  "transcript_is_canonical" %in% biomaRt::listAttributes(mart = mart)$name
#Check if three_prime_activated is activated or else if transcipts are cononical
if (has_canonical && three_prime_activated) {
  attributes <- c(attributes, "transcript_is_canonical", "chromosome_name",
    "transcript_mane_select", "ensembl_transcript_id_version")
  has_mane_select <-
  "transcript_mane_select" %in% biomaRt::listAttributes(mart = mart)$name
}else if (has_canonical) {
     attributes <- c(attributes, "transcript_is_canonical")
}
t2g <- biomaRt::getBM(
attributes = attributes,
mart = mart,
useCache = FALSE
)
# Set columns as NA if three_prime_activated is set to false or if the transcipts are not canonical
if (!has_canonical || !three_prime_activated) {
  t2g <- t2g %>% add_column(chromosome_name = NA, transcript_mane_select = NA,
      ensembl_transcript_id_version = NA)
}else if (!has_canonical) {
   t2g <- t2g %>% add_column(transcript_is_canonical = NA)
}
t2g <- t2g %>%
  rename(
    target_id = ensembl_transcript_id,
    ens_gene = ensembl_gene_id,
    ext_gene = external_gene_name,
    gene_desc = description,
    canonical = transcript_is_canonical,
    chromosome_name = chromosome_name,
    transcript_mane_select = transcript_mane_select,
    ensembl_transcript_id_version = ensembl_transcript_id_version,
  ) %>%
  mutate_at(
    vars(gene_desc),
    function(values) {
      str_trim(map(values, function(v) {
        str_split(v, r"{\[}")[[1]][1]
      }))
    } # remove trailing source annotation (e.g. [Source:HGNC Symbol;Acc:HGNC:5])
  ) %>%
  mutate_at(
    vars(canonical),
    function(values) {
      as_vector(
        map(
          str_trim(values),
          function(v) {
            if (is.na(v)) {
              NA
            } else if (v == "1") {
              TRUE
            } else if (v == "0") {
              FALSE
            }
          }
        )
      )
    }
  )
# Check if 3-prime-rna-seq is activated, filter transcipts that are mane selected and filter chromosomes that are defined as "patch" 
if (three_prime_activated) {
  if (has_mane_select) {
    t2g <- t2g %>%
    filter(!str_detect(chromosome_name, "patch|PATCH")) %>%
    filter(str_detect(transcript_mane_select, ""))
  }else {
    stop(
      str_c(
        "needed mane_selected column in biomart if three prime mode is activated"
      )
    )
  }
}
write_rds(t2g, file = snakemake@output[[1]], compress = "gz")