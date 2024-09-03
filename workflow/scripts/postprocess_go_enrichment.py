# import pandas as pd
from ast import literal_eval as make_tuple
import polars as pl


def calculate_enrichment(ratio_str):
    """Calculate enrichment based on ratios provided as strings."""
    ratio = make_tuple(ratio_str)
    enrichment = ratio[0] / ratio[1]
    return enrichment


def extract_study_items(value):
    if value and value.strip() != "":
        gene_values = value.split(", ")
        data = []
        for gene_value in gene_values:
            parts = gene_value.split(":")
            gene = parts[0]
            val = float(parts[1])
            data.append({"gene": gene, "value": val})
        return data
    else:
        return []


def calculate_sums(parsed_terms):
    return sum(abs(item["value"]) for item in parsed_terms)


# Load data
df_enr = pl.read_csv(snakemake.input["enrichment"], separator="\t")
df_sig = pl.read_csv(snakemake.input["significant_terms"], separator="\t")

# We only want to have the study items from df_sig
df_sig = df_sig.select(["GO", "study_items"])
df_merged = df_enr.join(df_sig, on="GO", how="inner", suffix="_df2")
df_merged = df_merged.with_columns(pl.col("study_items_df2").alias("study_items")).drop(
    "study_items_df2"
)


if not df_merged.is_empty():
    df_merged = df_merged.with_columns(
        pl.struct(
            [
                pl.col("ratio_in_study")
                .map_elements(calculate_enrichment)
                .alias("enrichment_study"),
                pl.col("ratio_in_pop")
                .map_elements(calculate_enrichment)
                .alias("enrichment_pop"),
            ]
        )
        .map_elements(
            lambda row: f"({row['enrichment_study']}, {row['enrichment_pop']})"
        )
        .alias("enrichment")
    )
else:
    df_merged = df_merged.with_columns(enrichment=pl.lit(None))

# The effect size is the sum of absolute values of the study items
df = df_merged.with_columns(
    [
        pl.col("study_items")
        .map_elements(lambda x: calculate_sums(extract_study_items(x)))
        .alias("effect")
    ]
)


df = df.with_columns(
    (-pl.col("p_fdr_bh").log(base=10) * pl.col("effect")).alias("signed_pi_value")
)

df_sorted = df.sort(pl.col("signed_pi_value").abs(), descending=True)

# Save the result to a file
df_sorted.write_csv(snakemake.output[0], separator="\t")
