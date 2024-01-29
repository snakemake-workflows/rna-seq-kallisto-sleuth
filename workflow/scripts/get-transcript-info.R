 log <- file(snakemake@log[[1]], open="wt")
 sink(log)
 sink(log, type="message")

library("biomaRt")
# tidy syntax
library("tidyverse")
# useful error messages upon aborting
library("cli")

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
        cli_abort(
          str_c(
            "Have tried all 4 available Ensembl biomaRt mirrors ",
            rounds,
            " times. You might have a connection problem, or no mirror is responsive.\n",
            "The last error message was:\n",
            message(e)
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

# define and keep those separately, to filter out below
sleuth_attributes <- c(
  "ensembl_transcript_id",
  "ensembl_gene_id",
  "external_gene_name",
  "description"
)

wanted_attributes <- sleuth_attributes

# get attributes to be able to check below, whether this species and version has
# ensembl canonical and MANE select annotations available
available_attributes <- biomaRt::listAttributes(mart = mart)$name

use_if_available <- function(attribute_name, available_attributes) {
  if (attribute_name %in% available_attributes) {
    attribute_name
  }
}

wanted_attributes <- c(
  wanted_attributes,
  use_if_available("transcript_is_canonical", available_attributes),
  # always get this if present, as this might be a useful annotation
  # for final results
  use_if_available("transcript_mane_select", available_attributes)
)


three_prime_attributes <- c(
  "ensembl_transcript_id_version",
  "chromosome_name",
  "transcript_length",
  "strand"
)

if (three_prime_activated & 
  !("transcript_mane_select" %in% available_attributes) &
  !("transcript_is_canonical" %in% available_attributes)
  ) {
  cli_abort(
    str_c(
      "Three prime mode for Lexogen QuantSeq analysis is activated, which ",
      "needs the transcript_mane_select or the transcript_is_canonical ",
      "attribute from biomart. However, these attributes ",
      "are not available for the species '",
      snakemake@params[["species"]],
      "' in the ensembl release version: ",
      snakemake@params[["version"]]
    )
  )
}

wanted_attributes <- c(
  wanted_attributes,
  three_prime_attributes
)

all_annotations <- biomaRt::getBM(
  attributes = wanted_attributes,
  mart = mart,
  useCache = FALSE
) |> as_tibble()


column_renames <- c(
  target_id = "ensembl_transcript_id",
  ens_gene = "ensembl_gene_id",
  ext_gene = "external_gene_name",
  gene_desc = "description",
  canonical = "transcript_is_canonical",
  mane = "transcript_mane_select"
)

t2g <- all_annotations |>
  rename(
    any_of(
      column_renames
    )
  ) |>
  select(
    -any_of(three_prime_attributes)
  ) |>
  mutate(
    # remove trailing source annotation (e.g. " [Source:HGNC Symbol;Acc:HGNC:5]")
    gene_desc = str_replace(
      gene_desc,
      " +\\[[^\\[\\]]+\\]",
      ""
    ),
    across(
      any_of("canonical"),
      ~ case_match(
        .x,
        NA ~ NA,
        1 ~ TRUE,
        0 ~ FALSE
      )
    )
  )

write_rds(
  t2g,
  file = snakemake@output[[1]],
  compress = "gz"
)

other_annotations <- all_annotations |>
  # TODO: determine why this filtering is done, and if this should also happen
  # for non-3-prime kallisto-sleuth input
  filter(!str_detect(chromosome_name, "patch|PATCH")) |>
  select(-c(chromosome_name, sleuth_attributes)) |>
  rename(transcript = ensembl_transcript_id_version) |>
  mutate(
    strand = case_match(
        strand,
        1 ~ "+",
        -1 ~ "-"
    )
  )

if ("transcript_mane_select" %in% colnames(other_annotations)) {
  other_annotations <- other_annotations |>
    select(
      -any_of("transcript_is_canonical")
    ) |>
    mutate(
      # we don't need the NCBI IDs that match the MANE transcripts, only an
      # indicator whether a transcript is in MANE select
      main_transcript_per_gene = if_else(transcript_mane_select != "", 1, 0, NA)
    ) |>
    select(
      -any_of("transcript_mane_select")
    )
} else {
  other_annotations <- other_annotations |>
    rename(
    # ensure consistent column presence in the output file
      main_transcript_per_gene = "transcript_is_canonical"
    )
}

write_tsv(
  other_annotations,
  snakemake@output[[2]]
)
