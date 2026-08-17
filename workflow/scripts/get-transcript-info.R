 log <- file(snakemake@log[[1]], open="wt")
 sink(log)
 sink(log, type="message")

library("biomaRt")
# tidy syntax
library("tidyverse")
# useful error messages upon aborting
library("cli")

# As ensembl has undergone extensive changes the get-transcript-info
# is now set to the look up the stable archive URL for the requested version
# Static lookup, because of unstable site during Ensembl's platform migration in 2026
archive_hosts <- c(
  "116" = "https://jun2026.archive.ensembl.org",
  "115" = "https://sep2025.archive.ensembl.org",
  "114" = "https://may2025.archive.ensembl.org",
  "113" = "https://oct2024.archive.ensembl.org",
  "112" = "https://may2024.archive.ensembl.org",
  "111" = "https://jan2024.archive.ensembl.org",
  "110" = "https://jul2023.archive.ensembl.org",
  "109" = "https://feb2023.archive.ensembl.org",
  "108" = "https://oct2022.archive.ensembl.org",
  "107" = "https://jul2022.archive.ensembl.org"
)

requested_version <- as.character(snakemake@params[["version"]])
host_url <- archive_hosts[[requested_version]]

if (is.null(host_url)) {
  cli_abort(str_c(
    "No archive host configured for Ensembl version ", requested_version,
    ". Known versions: ", str_c(names(archive_hosts), collapse = ", ")
  ))
}

rounds <- 0
mart <- NULL
while (is.null(mart) || class(mart)[[1]] != "Mart") {
  mart <- tryCatch(
    {
      rounds <- rounds + 1
      biomaRt::useEnsembl(
        biomart = "ENSEMBL_MART_ENSEMBL",
        dataset = str_c(snakemake@params[["species"]], "_gene_ensembl"),
        host    = host_url
      )
    },
    error = function(e) {
      if (rounds >= 3) {
        cli_abort(str_c(
          "Failed to connect to the Ensembl ", requested_version,
          " archive (", host_url, ") after ", rounds, " tries. Last error:\n",
          conditionMessage(e)
        ))
      }
      Sys.sleep(30)
      NULL
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
  use_if_available("ensembl_transcript_id_version", available_attributes),
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
write_tsv(
  t2g,
  file = snakemake@output[[2]],
)

other_annotations <- all_annotations |>
  # TODO: determine why this filtering is done, and if this should also happen
  # for non-3-prime kallisto-sleuth input
  filter(!str_detect(chromosome_name, "patch|PATCH")) |>
  select(-c(chromosome_name, sleuth_attributes)) |>
  rename(any_of(c(transcript = "ensembl_transcript_id_version"))) |>
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
  snakemake@output[[3]]
)
