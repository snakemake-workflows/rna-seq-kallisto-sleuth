import pandas as pd


# Read the TSV files
logcount_matrix = pd.read_csv(snakemake.input['logcount'], sep='\t')
genes_representative = pd.read_csv(
    snakemake.input['genes_representative'], sep='\t')

# Filter logcount_matrix to only include rows where 'transcript' is in 'target_id' of genes_representative
filtered_logcount_matrix = logcount_matrix[logcount_matrix['transcript'].isin(
    genes_representative['target_id'])]

# Save the filtered dataframe to a new TSV file
filtered_logcount_matrix.to_csv(snakemake.output[0], sep='\t', index=False)
