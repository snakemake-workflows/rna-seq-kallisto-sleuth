import sys
import pandas as pd

sys.stderr = open(snakemake.log[0], "w")

# Load input files
tpm_df = pd.read_csv(snakemake.input["tpm"], sep="\t")
diffexp_df = pd.read_csv(snakemake.input["diffexp"], sep="\t")

# Match on 'transcript' column
tpm_df["__sort_index"] = tpm_df["transcript"].map(
    {transcript: i for i, transcript in enumerate(diffexp_df["target_id"])}
)

# Unmatched transcripts get NaN -> pushed to end
tpm_df = tpm_df.sort_values(
    by="__sort_index", na_position="last"
).drop(columns="__sort_index")

# Write output
tpm_df.to_csv(snakemake.output[0], sep="\t", index=False)
