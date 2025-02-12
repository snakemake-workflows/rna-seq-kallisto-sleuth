import pandas as pd
import sys

sys.stderr = open(snakemake.log[0], "w")

# Read the TSV files
logcount_matrix = pd.read_csv(snakemake.input["logcount"], sep="\t")
diffexp = pd.read_csv(snakemake.input["diffexp"], sep="\t")

# Filter logcount_matrix to only include rows where 'transcript' is in 'target_id' of genes_representative
filtered_logcount_matrix = logcount_matrix[
    logcount_matrix["transcript"].isin(diffexp["target_id"])
]

# Save the filtered dataframe to a new TSV file
filtered_logcount_matrix.to_csv(snakemake.output[0], sep="\t", index=False)
