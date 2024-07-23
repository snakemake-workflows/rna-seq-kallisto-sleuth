import pandas as pd

# Load data
df = pd.read_csv(snakemake.input["spia"], sep='\t')

# Sortieren des DataFrames nach der Spalte 'total perturbation accumulation'
df_sorted = df.sort_values(
    by='total perturbation accumulation', ascending=True)

# Save the result to a file
df_sorted.to_csv(snakemake.output[0], sep='\t', index=False)
