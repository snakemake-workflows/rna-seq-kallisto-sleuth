import pandas as pd

df = pd.read_csv(snakemake.input["spia_table"], sep='\t')
df['gene_ratio'] = "(" + df['number of DE genes per pathway'].astype(str) + \
    ', ' + df['number of genes on the pathway'].astype(str) + ")"
df_sorted = df.sort_values(by='Combined FDR')
df_sorted.to_csv(snakemake.output[0], sep='\t', index=False)
