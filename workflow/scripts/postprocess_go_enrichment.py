# import pandas as pd
from ast import literal_eval as make_tuple
import polars as pl



def calculate_enrichment(ratio_str):
    """Calculate enrichment based on ratios provided as strings."""
    ratio = make_tuple(ratio_str)
    enrichment = ratio[0] / ratio[1]
    return enrichment


def sort_group(group):
    """Sort a DataFrame group by the 'p_uncorrected' column in ascending order."""
    return group.sort_values(by='p_uncorrected', ascending=True)


# Load data
df_enr = pl.read_csv(snakemake.input["enrichment"], separator='\t')
df_sig = pl.read_csv(snakemake.input["significant_terms"], separator='\t')

# We only want to have the study items from df_sig
df_sig = df_sig.select(["GO", "study_items"])
df_merged = df_enr.join(df_sig, on="GO", how="inner", suffix="_df2")
df_merged_sorted = df_merged.with_columns(
    pl.col("study_items_df2").alias("study_items")
).drop("study_items_df2")


if not df_merged_sorted.is_empty():
    df_merged_sorted = df_merged_sorted.with_columns(
        pl.struct([
            pl.col("ratio_in_study").map_elements(calculate_enrichment).alias("enrichment_study"),
            pl.col("ratio_in_pop").map_elements(calculate_enrichment).alias("enrichment_pop")
        ]).map_elements(lambda row: f"({row['enrichment_study']}, {row['enrichment_pop']})").alias("enrichment")
    )
else:
    df_merged_sorted = df_merged_sorted.with_columns(
        enrichment = pl.lit(None)
    )

# Save the result to a file
df_merged_sorted.write_csv(snakemake.output[0], separator='\t')
